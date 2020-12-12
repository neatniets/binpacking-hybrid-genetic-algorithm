#include "bin-packing.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <string.h>

static bool USE_CASE_INJECTION = false;
static enum init_type INIT = SUCCESSIVE_MUT;
static bool USE_LOCAL_SEARCH = true;
static enum search_type SEARCH = SWAP_RAND;
static enum search_adaptation_type ADAPT = LAMARCKIAN;
static int MAX_GENERATIONS
        // = INT_MAX;
        = 1000;
static double MAX_TIME
        = DBL_MAX;
        // = 3;

static const int MAX_THREADS = 4;

#define _STR(TOK)       #TOK
#define STR(TOK)        _STR(TOK)
#define CHECK(NUM, VAL1, VAL2, TARGET) \
        do { \
                if (strcmp(argv[NUM], STR(VAL1)) == 0) { \
                        TARGET = VAL1; \
                } else if (strcmp(argv[NUM], STR(VAL2)) == 0) { \
                        TARGET = VAL2; \
                } else { \
                        fprintf(stderr, STR(VAL1) " or " STR(VAL2) \
                                        " for argument " STR(NUM) "\n"); \
                        return -1; \
                } \
        } while(0)

int main(int argc, char **argv) {
        if (argc != 7) {
                fprintf(stderr, "bad number of args: %d\n", argc);
                return -1;
        }
        CHECK(1, true, false, USE_CASE_INJECTION);
        CHECK(2, SUCCESSIVE_MUT, HILL_CLIMB, INIT);
        CHECK(3, true, false, USE_LOCAL_SEARCH);
        CHECK(4, SWAP_RAND, SHUFFLE_GROUPS, SEARCH);
        CHECK(5, LAMARCKIAN, BALDWINIAN, ADAPT);
        if (argv[6][0] == '1') {
                MAX_GENERATIONS = INT_MAX;
                MAX_TIME = 3;
        } else if (argv[6][0] == '0') {
                MAX_GENERATIONS = 1000;
                MAX_TIME = DBL_MAX;
        } else {
                fprintf(stderr, "1 or 0 for argument 6\n");
                return -1;
        }
        size_t num_problems;
        scanf(" %zu", &num_problems);
        for (size_t i=0; i<num_problems; i++) {
                printf("PROBLEM #%zu:\n", i);
                /* skip problem identifier */
                scanf(" %*s");
                double bin_cap;
                size_t inst_sz, optimal_num_bins;
                scanf(" %lf %zu %zu",
                      &bin_cap, &inst_sz, &optimal_num_bins);
                double *prob_inst = malloc(inst_sz * sizeof(*prob_inst));
                if (prob_inst == NULL) {
                        abort();
                }
                for (size_t i=0; i<inst_sz; i++) {
                        scanf(" %lf", prob_inst+i);
                }

                struct solution sol;
                sol = genetic_algorithm(prob_inst,
                                        inst_sz,
                                        bin_cap,
                                        USE_CASE_INJECTION,
                                        INIT,
                                        USE_LOCAL_SEARCH,
                                        SEARCH,
                                        ADAPT,
                                        MAX_THREADS,
                                        MAX_GENERATIONS,
                                        MAX_TIME,
                                        stdout);
                solution_destroy(sol);
                free(prob_inst);
        }
        return 0;
}

static bool is_char_in_str(char ch, const char *str) {
        while (*str != '\0') {
                if (ch == *str) {
                        return true;
                }
                str++;
        }
        return false;
}
