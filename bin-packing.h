#ifndef BIN_PACKING_H
#define BIN_PACKING_H

#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include "bp-solution.h"

enum init_type {
        SUCCESSIVE_MUT,
        HILL_CLIMB
};
enum search_type {
        NONE,
        SWAP_RAND,
        SHUFFLE_GROUPS,
        DOMINANCE
};
enum search_adaptation_type {
        LAMARCKIAN,
        BALDWINIAN
};

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
                                  FILE *out);

#endif /* !BIN_PACKING_H */
