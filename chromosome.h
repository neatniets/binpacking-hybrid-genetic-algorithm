#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "bp-solution.h"
#include <stdbool.h>

struct chromosome {
        double fitness;
        size_t *perm;
        struct solution sol;
};
struct bald_chrom {
        struct chromosome chrom;
        struct solution bald_sol;
};

void chrom_init(struct chromosome *chrom,
                bool is_baldwinian);
void chrom_destroy(struct chromosome *chrom,
                   bool is_baldwinian);
void chrom_eval(struct chromosome *chrom,
                const double *prob_inst,
                size_t inst_sz,
                double bin_cap);

void chrom_mut(struct chromosome *chrom,
               size_t inst_sz);
struct chromosome chrom_cx(struct chromosome parent1,
                           struct chromosome parent2,
                           size_t inst_sz);

struct search_flags {
        bool perm_modified : 1;
        bool sol_modified : 1;
};
typedef struct search_flags (*chrom_search_func)(size_t *perm,
                                                 struct solution *sol,
                                                 const double *prob_inst,
                                                 size_t inst_sz,
                                                 double bin_cap);

/* returns number of searches conducted */
int chrom_search(struct chromosome *chrom,
                 bool is_baldwinian,
                 const double *prob_inst,
                 size_t inst_sz,
                 double bin_cap,
                 bool is_greedy,
                 int max_searches,
                 chrom_search_func get_neighbor);

struct search_flags chrom_search_swap(size_t *perm,
                                      struct solution *unused,
                                      const double *prob_inst,
                                      size_t inst_sz,
                                      double bin_cap);
struct search_flags chrom_search_shuffle(size_t *unused,
                                         struct solution *sol,
                                         const double *prob_inst,
                                         size_t inst_sz,
                                         double bin_cap);
struct search_flags chrom_search_dom(size_t *unused,
                                     struct solution *sol,
                                     const double *prob_inst,
                                     size_t inst_sz,
                                     double bin_cap);

#endif /* !CHROMOSOME_H */
