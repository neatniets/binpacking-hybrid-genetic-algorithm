#ifndef BP_SOLUTION_H
#define BP_SOLUTION_H

#include <stddef.h>
#include <stdio.h>

struct bin {
        double item_sum;
        int num_items;
        size_t *item_indices;
};
void bin_init(struct bin *restrict bin);
void bin_add(struct bin *restrict bin,
             size_t item_index,
             const double *restrict prob_inst);
void bin_copy(struct bin *restrict dest,
              const struct bin src);
void bin_print(struct bin bin,
               FILE *restrict out);

struct solution {
        int num_bins;
        struct bin *bins;
};
void solution_init(struct solution *restrict sol);
void solution_add(struct solution *restrict sol,
                  struct bin bin);
void solution_destroy(struct solution sol);
void solution_copy(struct solution *restrict dest,
                   const struct solution src);

void solution_first_fit(struct solution *restrict sol,
                        const double *restrict prob_inst,
                        size_t inst_sz,
                        const size_t *restrict perm,
                        double bin_cap);
size_t *solution_reverse_first_fit(struct solution sol,
                                   size_t perm_sz);
void solution_print(struct solution sol,
                    FILE *restrict out);

double solution_eval(struct solution sol,
                     double bin_cap);

#endif /* !BP_SOLUTION_H */
