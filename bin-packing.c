#include "bin-packing.h"
#include "chromosome.h"
#include <assert.h>
#include <stdlib.h>
#include "parallel-foreach.h"
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

static FILE *CASE_INJECT_FILE = NULL;

#define POP_I(POP, I) \
        (*(struct chromosome *)((POP).chroms \
                                + ((I) * (((POP).is_baldwinian) \
                                          ? sizeof(struct bald_chrom) \
                                          : sizeof(struct chromosome)))))
struct population {
        void *chroms;
        bool is_baldwinian;
        size_t pop_sz;
        double bin_cap;
        const double *prob_inst;
        size_t inst_sz;
};

static double sum_inst(const double *prob_inst, size_t inst_sz);
static double time_elapsed(const struct timespec *time_start);
struct foreach_init_context {
        const bool is_baldwinian;
        const size_t perm_sz;
};
static int pop_init_foreach_init(void *elem,
                                 void *context);
static void pop_init_mut(struct chromosome *chrom,
                         size_t inst_sz);
/* returns number of searches conducted, if any */
static int pop_init(struct population pop,
                    enum init_type init,
                    chrom_search_func search_func,
                    bool use_case_injection,
                    int max_threads);
struct eval_foreach_context {
        const struct population pop;
};
static int pop_eval_foreach(void *elem,
                            void *eval_foreach_context);
static void pop_eval(struct population pop,
                     int max_threads);
struct select_foreach_context {
        const struct population pop;
};
static int tournament_select_foreach(void *elem,
                                     void *select_foreach_context);
struct cx_foreach_context {
        const struct population pop;
        const size_t *tourn;
        const size_t tourn_count;
};
static int pop_cx_foreach(void *elem,
                          void *cx_foreach_context);
static void pop_cx(struct population pop,
                   struct population *children,
                   int max_threads);
static void pop_replace(struct population pop,
                        struct population children,
                        int max_threads);
static void destroy_children(struct population *children,
                             int max_threads);
struct mut_foreach_context {
        const double mut_rate;
        const size_t inst_sz;
};
static int pop_mut_foreach(void *elem,
                           void *mut_foreach_context);
static void pop_mut(struct population pop,
                    double mut_rate,
                    int max_threads);
struct search_foreach_context {
        const struct population pop;
        const chrom_search_func search_func;
        int num_searches;
        volatile bool is_locked;
};
static int pop_search_foreach(void *elem,
                              void *search_foreach_context);
/* returns number of searches conducted, if any */
static int pop_search(struct population pop,
                      chrom_search_func search_func,
                      int max_threads);
struct avg_fitness_foreach_context {
        double sum_fitness;
        volatile bool is_locked;
};
static int pop_avg_fitness_foreach(void *elem,
                                   void *avg_fitness_foreach_context);
static double pop_avg_fitness(const struct population pop,
                              int max_threads);
struct fitness_foreach_context {
        double best_fitness;
        struct chromosome *best_chrom;
        volatile bool is_locked;
};
static int pop_best_fitness_foreach(void *elem,
                                    void *fitness_foreach_context);
static double pop_best_fitness(const struct population pop,
                               struct solution *best_sol_copy,
                               int max_threads);

struct solution genetic_algorithm(const double *prob_inst,
                                  size_t inst_sz,
                                  double bin_cap,
                                  bool use_case_injection,
                                  enum init_type init,
                                  bool use_local_search,
                                  enum search_type search,
                                  enum search_adaptation_type adapt,
                                  int max_threads,
                                  int max_generations,
                                  double max_time,
                                  FILE *out) {
        if (!use_local_search && (adapt == BALDWINIAN)) {
                assert(false);
        }
        const int pop_sz = 100;
        const double mut_rate = 0.1;
        struct timespec time_start;
        clock_gettime(CLOCK_REALTIME, &time_start);
        const int theoretical_min_bins = sum_inst(prob_inst, inst_sz)
                                         / bin_cap;
        int num_searches;
        int num_gens = 1;
        double best_fitness;
        double avg_fitness;
        struct solution best_sol;

        chrom_search_func search_func;
        switch (search) {
                case NONE:
                        search_func = NULL;
                        break;
                case SWAP_RAND:
                        search_func = chrom_search_swap;
                        break;
                case SHUFFLE_GROUPS:
                        search_func = chrom_search_shuffle;
                        break;
                case DOMINANCE:
                        search_func = chrom_search_dom;
                        break;
                default:
                        assert(false);
        }

        struct population pop;
        pop.pop_sz = pop_sz;
        pop.bin_cap = bin_cap;
        pop.prob_inst = prob_inst;
        pop.inst_sz = inst_sz;
        switch (adapt) {
                case LAMARCKIAN:
                        pop.chroms = malloc(pop_sz
                                            * sizeof(struct chromosome));
                        pop.is_baldwinian = false;
                        break;
                case BALDWINIAN:
                        pop.chroms = malloc(pop_sz
                                            * sizeof(struct bald_chrom));
                        pop.is_baldwinian = true;
                        break;
                default:
                        assert(false);
        }
        if (pop.chroms == NULL) {
                abort();
        }

        num_searches = pop_init(pop, init, search_func,
                                use_case_injection, max_threads);
        pop_eval(pop, max_threads);
        solution_init(&best_sol);
        best_fitness = pop_best_fitness(pop, &best_sol, max_threads);
        avg_fitness = pop_avg_fitness(pop, max_threads);

        fprintf(out, "theoretical minimum bins: %d\n"
               "gen\tcum time\tbes bin\tbest fit\tavrg fit\tsearches\n",
               theoretical_min_bins);
        fprintf(out, "%d\t%lf\t%d\t%lf\t%lf\t%d\n",
               num_gens, time_elapsed(&time_start), best_sol.num_bins,
               best_fitness, avg_fitness, num_searches);

        num_gens++;

        while ((num_gens <= max_generations)
               && (time_elapsed(&time_start) <= max_time)
               && (best_sol.num_bins > theoretical_min_bins)) {
                struct population children;
                pop_cx(pop, &children, max_threads);
                pop_replace(pop, children, max_threads);
                destroy_children(&children, max_threads);

                if (!use_local_search) {
                        pop_mut(pop, mut_rate, max_threads);
                } else {
                        num_searches += pop_search(pop, search_func,
                                                   max_threads);
                }

                pop_eval(pop, max_threads);
                solution_destroy(best_sol);
                best_fitness = pop_best_fitness(pop, &best_sol, max_threads);
                avg_fitness = pop_avg_fitness(pop, max_threads);
                fprintf(out, "%d\t%lf\t%d\t%lf\t%lf\t%d\n",
                       num_gens, time_elapsed(&time_start), best_sol.num_bins,
                       best_fitness, avg_fitness, num_searches);

                num_gens++;
        }

        for (size_t i = 0; i < pop_sz; i++) {
                chrom_destroy(&POP_I(pop, i), pop.is_baldwinian);
        }
        free(pop.chroms);
        if (use_case_injection) {
                size_t *buf = solution_reverse_first_fit(best_sol, inst_sz);
                fseek(CASE_INJECT_FILE, 0, SEEK_END);
                fwrite(buf, sizeof(*buf), inst_sz, CASE_INJECT_FILE);
                free(buf);
        }
        return best_sol;
}

static double sum_inst(const double *prob_inst, size_t inst_sz) {
        double sum = 0.0;
        for (size_t i = 0; i < inst_sz; i++) {
                sum += prob_inst[i];
        }
        return sum;
}
static double time_elapsed(const struct timespec *time_start) {
        struct timespec new_time;
        clock_gettime(CLOCK_REALTIME, &new_time);
        double sec = new_time.tv_sec - time_start->tv_sec;
        double nsec = (long)new_time.tv_nsec - time_start->tv_nsec;
        nsec /= 1000000000;
        return sec + nsec;
}
static int pop_init_foreach_init(void *elem,
                                 void *context) {
        struct foreach_init_context *con = context;
        struct chromosome *chrom = elem;
        chrom_init(chrom, con->is_baldwinian);
        chrom->perm = malloc(con->perm_sz * sizeof(*chrom->perm));
        if (chrom->perm == NULL) {
                abort();
        }
        return 0;
}
static void pop_init_mut(struct chromosome *chrom,
                         size_t inst_sz) {
        const int max_muts = 20;
        const int num_muts = (rand() % max_muts) + 1;
        for (int i = 0; i < num_muts; i++) {
                chrom_mut(chrom, inst_sz);
        }
}
/* returns number of searches conducted, if any */
static int pop_init(struct population pop,
                    enum init_type init,
                    chrom_search_func search_func,
                    bool use_case_injection,
                    int max_threads) {
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? sizeof(struct bald_chrom)
                                  : sizeof(struct chromosome);
        /* initialize all the chromosomes */
        struct foreach_init_context tmp1 = {.is_baldwinian = pop.is_baldwinian,
                                            .perm_sz = pop.inst_sz};
        parallel_foreach(max_threads, pop.chroms, pop.pop_sz, chrom_size,
                         &tmp1, pop_init_foreach_init);

        /* read saved data from case injection; THE PROGRAM IS ASSUMED TO
         * ONLY RUN ON DATA OF EQUAL LENGTHS UNTIL TERMINATION */
        size_t i = 0;
        if (use_case_injection) {
                /* initialize temporary file once for entire run of program */
                if (CASE_INJECT_FILE == NULL) {
                        CASE_INJECT_FILE = tmpfile();
                        if (CASE_INJECT_FILE == NULL) {
                                abort();
                        }
                }
                /* read data from temporary file as it was inserted */
                rewind(CASE_INJECT_FILE);
                while (i < pop.pop_sz) {
                        size_t bytes_read = fread(POP_I(pop, i).perm,
                                                  sizeof(size_t), pop.inst_sz,
                                                  CASE_INJECT_FILE);
                        if (bytes_read == 0) {
                                if (feof(CASE_INJECT_FILE)) {
                                        break;
                                } else if (ferror(CASE_INJECT_FILE)) {
                                        abort();
                                }
                        }
                        i++;
                }
        }

        /* if there is no previous result, make one */
        if (i == 0) {
                size_t *perm = POP_I(pop, i).perm;
                for (size_t j = 0; j < pop.inst_sz; j++) {
                        perm[j] = j;
                }
                pop_init_mut(&POP_I(pop, i), pop.inst_sz);
                i++;
        }
        if (init == SUCCESSIVE_MUT) { // mutate the previous result
                for (; i < pop.pop_sz; i++) {
                        memcpy(POP_I(pop, i).perm, POP_I(pop, i-1).perm,
                               pop.inst_sz * sizeof(*POP_I(pop, i).perm));
                        pop_init_mut(&POP_I(pop, i), pop.inst_sz);
                }
                return 0;
        } else if (init == HILL_CLIMB) { // use hill-climbing to initialize
                const int max_searches = 100;
                int num_searches = 0;
                for (; i < pop.pop_sz; i++) {
                        memcpy(POP_I(pop, i).perm, POP_I(pop, i-1).perm,
                               sizeof(*POP_I(pop, i).perm) * pop.inst_sz);
                        chrom_eval(&POP_I(pop, i), pop.prob_inst,
                                   pop.inst_sz, pop.bin_cap);
                        num_searches += chrom_search(&POP_I(pop, i), false,
                                                     pop.prob_inst,
                                                     pop.inst_sz, pop.bin_cap,
                                                     true, max_searches,
                                                     search_func);
                        if (memcmp(POP_I(pop, i).perm, POP_I(pop, i-1).perm,
                                   sizeof(*POP_I(pop, i).perm) * pop.inst_sz)
                            == 0) {
                                pop_init_mut(&POP_I(pop, i), pop.inst_sz);
                        }
                }
                return num_searches;
        } else {
                assert(false);
                return -1;
        }
}
static int pop_eval_foreach(void *elem,
                            void *eval_foreach_context) {
        struct eval_foreach_context *context = eval_foreach_context;
        struct chromosome *chrom = elem;
        chrom_eval(chrom, context->pop.prob_inst, context->pop.inst_sz,
                   context->pop.bin_cap);
        return 0;
}
static void pop_eval(struct population pop,
                     int max_threads) {
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? sizeof(struct bald_chrom)
                                  : sizeof(struct chromosome);
        parallel_foreach(max_threads, pop.chroms, pop.pop_sz,
                         chrom_size, &pop, pop_eval_foreach);
}
static int tournament_select_foreach(void *elem,
                                     void *select_foreach_context) {
        struct select_foreach_context *context = select_foreach_context;
        size_t *index = elem;
        size_t i1, i2;
        i1 = rand() % context->pop.pop_sz;
        while (i2 = rand() % context->pop.pop_sz, i2 == i1);
        if (POP_I(context->pop, i1).fitness
            > POP_I(context->pop, i2).fitness) {
                *index = i1;
        } else {
                *index = i2;
        }
        return 0;
}
static int pop_cx_foreach(void *elem,
                          void *cx_foreach_context) {
        struct cx_foreach_context *context = cx_foreach_context;
        struct chromosome *chrom = elem;
        size_t i1, i2;
        i1 = rand() % context->tourn_count;
        while (i2 = rand() % context->tourn_count, i2 == i1);
        *chrom = chrom_cx(POP_I(context->pop, i1), POP_I(context->pop, i2),
                          context->pop.inst_sz);
        return 0;
}
static void pop_cx(struct population pop,
                   struct population *children,
                   int max_threads) {
        pop_eval(pop, max_threads);

        /* tournament selection */
        const size_t tourn_count = pop.pop_sz;
        size_t *tourn = malloc(tourn_count * sizeof(*tourn));
        if (tourn == NULL) {
                abort();
        }
        parallel_foreach(max_threads, tourn, tourn_count, sizeof(*tourn),
                         &pop, tournament_select_foreach);

        /* initialize children */
        children->pop_sz = pop.pop_sz;
        children->is_baldwinian = false;
        children->chroms = malloc(children->pop_sz
                                  * sizeof(struct chromosome));
        struct chromosome *chroms = children->chroms;
        /* keep best chromosome */
        solution_init(&chroms[0].sol);
        pop_best_fitness(pop, &chroms[0].sol, max_threads);
        chroms[0].perm = solution_reverse_first_fit(chroms[0].sol,
                                                    pop.inst_sz);
        solution_destroy(chroms[0].sol);
        solution_init(&chroms[0].sol);
        chroms[0].fitness = -1;
        /* crossover */
        struct cx_foreach_context tmp1 = {.pop = pop,
                                          .tourn_count = tourn_count,
                                          .tourn = tourn};
        parallel_foreach(max_threads, chroms + 1, children->pop_sz - 1,
                         sizeof(*chroms), &tmp1, pop_cx_foreach);

        free(tourn);
}
static void pop_replace(struct population pop,
                        struct population children,
                        int max_threads) {
        for (size_t i = 0; i < pop.pop_sz; i++) {
                chrom_destroy(&POP_I(pop, i), pop.is_baldwinian);
                POP_I(pop, i) = POP_I(children, i);
        }
}
static void destroy_children(struct population *children,
                             int max_threads) {
        free(children->chroms);
}
static int pop_mut_foreach(void *elem,
                           void *mut_foreach_context) {
        struct mut_foreach_context *context = mut_foreach_context;
        struct chromosome *chrom = elem;
        const double roll = (double)rand() / RAND_MAX;
        if (roll <= context->mut_rate) {
                chrom_mut(chrom, context->inst_sz);
        }
        return 0;
}
static void pop_mut(struct population pop,
                    double mut_rate,
                    int max_threads) {
        struct mut_foreach_context tmp = {.mut_rate = mut_rate,
                                          .inst_sz = pop.inst_sz};
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? sizeof(struct bald_chrom)
                                  : sizeof(struct chromosome);
        parallel_foreach(max_threads, pop.chroms + (chrom_size),
                         pop.pop_sz - 1, chrom_size,
                         &tmp, pop_mut_foreach);
}
static int pop_search_foreach(void *elem,
                              void *search_foreach_context) {
        const int max_searches = 100;
        struct search_foreach_context *context = search_foreach_context;
        struct chromosome *chrom = elem;
        chrom_eval(chrom, context->pop.prob_inst, context->pop.inst_sz,
                   context->pop.bin_cap);
        int tmp = chrom_search(chrom, context->pop.is_baldwinian,
                               context->pop.prob_inst,
                               context->pop.inst_sz,
                               context->pop.bin_cap,
                               true, max_searches,
                               context->search_func);
        while (context->is_locked);
        context->is_locked = true;
        context->num_searches += tmp;
        context->is_locked = false;
        return 0;
}
/* returns number of searches conducted, if any */
static int pop_search(struct population pop,
                      chrom_search_func search_func,
                      int max_threads) {
        if (search_func == NULL) {
                return 0;
        }
        struct search_foreach_context tmp = {.pop = pop,
                                             .is_locked = false,
                                             .num_searches = 0,
                                             .search_func = search_func};
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? (sizeof(struct bald_chrom))
                                  : (sizeof(struct chromosome));
        parallel_foreach(max_threads, pop.chroms, pop.pop_sz, chrom_size,
                         &tmp, pop_search_foreach);
        return tmp.num_searches;
}
static int pop_avg_fitness_foreach(void *elem,
                                   void *avg_fitness_foreach_context) {
        struct avg_fitness_foreach_context *context
                = avg_fitness_foreach_context;
        struct chromosome *chrom = elem;
        while (context->is_locked);
        context->is_locked = true;
        context->sum_fitness += chrom->fitness;
        context->is_locked = false;
        return 0;
}
static double pop_avg_fitness(const struct population pop,
                              int max_threads) {
        struct avg_fitness_foreach_context tmp = {.sum_fitness = 0,
                                                  .is_locked = 0};
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? sizeof(struct bald_chrom)
                                  : sizeof(struct chromosome);
        parallel_foreach(max_threads, pop.chroms, pop.pop_sz,
                         chrom_size, &tmp, pop_avg_fitness_foreach);
        return tmp.sum_fitness / pop.pop_sz;
}
static int pop_best_fitness_foreach(void *elem,
                                    void *fitness_foreach_context) {
        struct fitness_foreach_context *context = fitness_foreach_context;
        struct chromosome *chrom = elem;
        while (context->is_locked); // mutex lock
        context->is_locked = true;
        if (chrom->fitness > context->best_fitness) {
                context->best_fitness = chrom->fitness;
                context->best_chrom = chrom;
        }
        context->is_locked = false;
        return 0;
}
static double pop_best_fitness(const struct population pop,
                               struct solution *best_sol_copy,
                               int max_threads) {
        struct fitness_foreach_context tmp = {.best_fitness = 0,
                                              .is_locked = false};
        const size_t chrom_size = (pop.is_baldwinian)
                                  ? sizeof(struct bald_chrom)
                                  : sizeof(struct chromosome);
        parallel_foreach(max_threads, pop.chroms, pop.pop_sz, chrom_size,
                         &tmp, pop_best_fitness_foreach);
        const struct solution *src_ptr;
        if (pop.is_baldwinian) {
                src_ptr = &((struct bald_chrom *)tmp.best_chrom)->bald_sol;
        } else {
                src_ptr = &tmp.best_chrom->sol;
        }
        solution_copy(best_sol_copy, *src_ptr);
        return tmp.best_fitness;
}
