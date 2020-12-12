#include "chromosome.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifndef NDEBUG
#include <stdio.h>
#endif

#define CHROM2BALD(CHROM_PTR) \
        ((struct bald_chrom *)CHROM_PTR)

static void perm_rand_swap(size_t *perm,
                           size_t perm_sz);
static size_t *perm_ox(const size_t *parent1,
                       const size_t *parent2,
                       size_t perm_sz);
static int sz_sort_asc(const void *a, const void *b);

void chrom_init(struct chromosome *chrom,
                bool is_baldwinian) {
        chrom->fitness = -1;
        chrom->perm = NULL;
        solution_init(&chrom->sol);
        if (is_baldwinian) {
                solution_init(&CHROM2BALD(chrom)->bald_sol);
        }
}
void chrom_destroy(struct chromosome *chrom,
                   bool is_baldwinian) {
        if (chrom == NULL) {
                return;
        }
        free(chrom->perm);
        solution_destroy(chrom->sol);
        if (is_baldwinian) {
                solution_destroy(CHROM2BALD(chrom)->bald_sol);
        }
}
void chrom_eval(struct chromosome *chrom,
                const double *prob_inst,
                size_t inst_sz,
                double bin_cap) {
        if (chrom->fitness >= 0) {
                return;
        }
        solution_destroy(chrom->sol);
        solution_first_fit(&chrom->sol, prob_inst, inst_sz,
                           chrom->perm, bin_cap);
        chrom->fitness = solution_eval(chrom->sol, bin_cap);
}

void chrom_mut(struct chromosome *chrom,
               size_t inst_sz) {
        if (chrom->fitness > 0) {
                chrom->fitness = -1;
                solution_destroy(chrom->sol);
                solution_init(&chrom->sol);
        }
        perm_rand_swap(chrom->perm, inst_sz);
}
/* OX - Order Crossover */
struct chromosome chrom_cx(struct chromosome parent1,
                           struct chromosome parent2,
                           size_t inst_sz) {
        struct chromosome child;
        chrom_init(&child, false);
        child.perm = perm_ox(parent1.perm, parent2.perm, inst_sz);
        return child;
}

int chrom_search(struct chromosome *chrom,
                 bool is_baldwinian,
                 const double *prob_inst,
                 size_t inst_sz,
                 double bin_cap,
                 bool is_greedy,
                 int max_searches,
                 chrom_search_func get_neighbor) {
        /* permutation to be passed to get_neighbor */
        size_t *working_perm = malloc(inst_sz * sizeof(*working_perm));
        if (working_perm == NULL) {
                abort();
        }
        /* solution struct to be passed to get_neighbor */
        struct solution working_sol;
        /* location to store current best solution */
        struct solution *best_sol_ptr;
        if (is_baldwinian) {
                best_sol_ptr = &CHROM2BALD(chrom)->bald_sol;
        } else {
                best_sol_ptr = &chrom->sol;
        }

        int num_searches = 0;
        for (; num_searches < max_searches; num_searches++) {
                /* revert changes that made a worse neighbor */
                memcpy(working_perm, chrom->perm,
                       inst_sz * sizeof(*working_perm));
                solution_copy(&working_sol, chrom->sol);

                /* flag returns used to determine if permutation needs to
                 * be recreated via reverse-first-fit or if solution needs
                 * to be created from permutation */
                struct search_flags flags;
                flags = get_neighbor(working_perm, &working_sol, prob_inst,
                                     inst_sz, bin_cap);
                if (!flags.sol_modified) {
                        solution_first_fit(&working_sol, prob_inst, inst_sz,
                                           working_perm, bin_cap);
                } else if (!flags.perm_modified) {
                        free(working_perm);
                        working_perm
                                = solution_reverse_first_fit(working_sol,
                                                             inst_sz);
                        solution_first_fit(&working_sol, prob_inst, inst_sz,
                                           working_perm, bin_cap);
                }

                double working_fitness = solution_eval(working_sol, bin_cap);
                if (working_fitness > chrom->fitness) {
                        /* accept new best permutation/solution */
                        chrom->fitness = working_fitness;
                        solution_destroy(*best_sol_ptr);
                        *best_sol_ptr = working_sol;
                        if (!is_baldwinian) {
                                memcpy(chrom->perm, working_perm,
                                       inst_sz * sizeof(*chrom->perm));
                        }
                        if (is_greedy) {
                                break;
                        }
                } else {
                        solution_destroy(working_sol);
                }
        }

        free(working_perm);
        return num_searches;
}

struct search_flags chrom_search_swap(size_t *perm,
                                      struct solution *unused,
                                      const double *prob_inst,
                                      size_t inst_sz,
                                      double bin_cap) {
        (void)unused;
        perm_rand_swap(perm, inst_sz);
        return (struct search_flags){.perm_modified = true,
                                     .sol_modified = false};
}
struct search_flags chrom_search_shuffle(size_t *unused,
                                         struct solution *sol,
                                         const double *prob_inst,
                                         size_t inst_sz,
                                         double bin_cap) {
        (void)unused;
        int i1, i2;
        i1 = rand() % sol->num_bins;
        while (i2 = rand() % sol->num_bins, i2 == i1);
        struct bin tmp = sol->bins[i1];
        sol->bins[i1] = sol->bins[i2];
        sol->bins[i2] = tmp;
        return (struct search_flags){.perm_modified = false,
                                     .sol_modified = true};
}
struct search_flags chrom_search_dom(size_t *unused,
                                     struct solution *sol,
                                     const double *prob_inst,
                                     size_t inst_sz,
                                     double bin_cap) {
        abort();
}

static void perm_rand_swap(size_t *perm,
                           size_t perm_sz) {
        size_t i1, i2;
        i1 = rand() % perm_sz;
        while (i2 = rand() % perm_sz, i2 == i1);
        size_t tmp = perm[i1];
        perm[i1] = perm[i2];
        perm[i2] = tmp;
}
static size_t *perm_ox(const size_t *parent1,
                       const size_t *parent2,
                       size_t perm_sz) {
        size_t *child = malloc(perm_sz * sizeof(*child));
        if (child == NULL) {
                abort();
        }

        /* cut points are just before i1/i2 */
        size_t i1, i2;
        i1 = rand() % (perm_sz + 1);
        while (i2 = rand() % (perm_sz + 1), i2 == i1);
        if (i1 > i2) {
                size_t tmp = i1;
                i1 = i2;
                i2 = tmp;
        }
        /* copy midsection of parent1 to child in same spot */
        memcpy(child + i1, parent1 + i1, (i2 - i1) * sizeof(*child));

        /* create sorted array of used values for binary search */
        size_t *used_vals = malloc((i2 - i1) * sizeof(*used_vals));
        if (used_vals == NULL) {
                abort();
        }
        memcpy(used_vals, parent1 + i1, (i2 - i1) * sizeof(*child));
        qsort(used_vals, (i2 - i1), sizeof(*used_vals), sz_sort_asc);

        /* copy all values not in child from parent2 starting at i2,
         * wrapping as necessary */
        size_t child_pos = i2;
        bool has_wrapped = false;
        for (size_t i = i2;
             (!(has_wrapped && (i >= i2)) && (child_pos != i1));
             i++) {
                if (i == perm_sz) {
                        i = 0;
                        has_wrapped = true;
                }
                /* search used_vals for current element */
                size_t *found_val = bsearch(parent2 + i, used_vals, (i2 - i1),
                                            sizeof(*found_val), sz_sort_asc);
                if (found_val == NULL) {
                        if (child_pos == perm_sz) {
                                child_pos = 0;
                        }
                        child[child_pos] = parent2[i];
                        child_pos++;
                }
        }
        free(used_vals);
        return child;
}
static int sz_sort_asc(const void *a, const void *b) {
        const size_t *av = a;
        const size_t *bv = b;
        if (*av > *bv) {
                return 1;
        } else if (*av < *bv) {
                return -1;
        } else {
                return 0;
        }
}
