#include "bp-solution.h"
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

static int bin_rff(struct bin bin,
                   size_t *restrict dest_perm);

void bin_init(struct bin *restrict bin) {
        *bin = (struct bin){.num_items = 0,
                            .item_indices = NULL,
                            .item_sum = 0.0};
}
void bin_add(struct bin *restrict bin,
             size_t item_index,
             const double *restrict prob_inst) {
        bin->num_items++;
        bin->item_indices = realloc(bin->item_indices,
                                    bin->num_items
                                    * sizeof(*bin->item_indices));
        bin->item_indices[bin->num_items - 1] = item_index;
        bin->item_sum += prob_inst[item_index];
}
void bin_copy(struct bin *restrict dest,
              const struct bin src) {
        *dest = (struct bin){.item_sum = src.item_sum,
                             .num_items = src.num_items,
                             .item_indices
                                     = malloc(src.num_items
                                              * sizeof(*dest->item_indices))};
        if (dest->item_indices == NULL) {
                abort();
        }
        memcpy(dest->item_indices, src.item_indices,
               dest->num_items * sizeof(*dest->item_indices));
}
void bin_print(struct bin bin,
               FILE *restrict out) {
        fprintf(out, "%lf | ", bin.item_sum);
        for (int i = 0; i < bin.num_items; i++) {
                fprintf(out, "%zu ", bin.item_indices[i]);
        }
}

void solution_init(struct solution *restrict sol) {
        *sol = (struct solution){.num_bins = 0,
                                 .bins = NULL};
}
void solution_add(struct solution *restrict sol,
                  struct bin bin) {
        sol->num_bins++;
        sol->bins = realloc(sol->bins,
                            sol->num_bins * sizeof(*sol->bins));
        sol->bins[sol->num_bins - 1] = bin;
}
void solution_destroy(struct solution s) {
        if (s.bins == NULL) {
                return;
        }
        for (int i = 0; i < s.num_bins; i++) {
                free(s.bins[i].item_indices);
        }
        free(s.bins);
}
void solution_copy(struct solution *restrict dest,
                   const struct solution src) {
        *dest = (struct solution){.num_bins = src.num_bins,
                                  .bins = malloc(src.num_bins
                                                 * sizeof(*dest->bins))};
        if (dest->bins == NULL) {
                abort();
        }
        for (int i = 0; i < dest->num_bins; i++) {
                bin_copy(dest->bins + i, src.bins[i]);
        }
}

void solution_first_fit(struct solution *restrict sol,
                        const double *restrict prob_inst,
                        size_t inst_sz,
                        const size_t *restrict perm,
                        double bin_cap) {
        if (sol->num_bins > 0) {
                solution_destroy(*sol);
                solution_init(sol);
        }
        for (size_t i = 0; i < inst_sz; i++) {
                int j = 0;
                for (; j < sol->num_bins; j++) {
                        if (sol->bins[j].item_sum + prob_inst[perm[i]]
                            <= bin_cap) {
                                bin_add(sol->bins + j, perm[i], prob_inst);
                                break;
                        }
                }
                if (j == sol->num_bins) {
                        struct bin tmp;
                        bin_init(&tmp);
                        bin_add(&tmp, perm[i], prob_inst);
                        solution_add(sol, tmp);
                }
        }
}
size_t *solution_reverse_first_fit(struct solution sol,
                                   size_t perm_sz) {
        size_t *perm = malloc(perm_sz * sizeof(*perm));
        size_t *cur_pos = perm;
        for (int i = 0; i < sol.num_bins; i++) {
                cur_pos += bin_rff(sol.bins[i], cur_pos);
        }
        return perm;
}
void solution_print(struct solution sol,
                    FILE *restrict out) {
        for (int i = 0; i < sol.num_bins; i++) {
                fprintf(out, "bin %d: ", i);
                bin_print(sol.bins[i], out);
                fputc('\n', out);
        }
}

double solution_eval(struct solution sol,
                     double bin_cap) {
        const double mul = 1.0 / (bin_cap * bin_cap * sol.num_bins);
        double sum = 0;
        for (int i = 0; i < sol.num_bins; i++) {
                sum += (sol.bins[i].item_sum * sol.bins[i].item_sum);
        }
        return sum * mul;
}

static int bin_rff(struct bin bin,
                   size_t *restrict dest_perm) {
        for (int i = 0; i < bin.num_items; i++) {
                dest_perm[i] = bin.item_indices[i];
        }
        return bin.num_items;
}
