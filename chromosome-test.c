#include "chromosome.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SEARCHES        100

static void chrom_print(const struct chromosome *chrom,
                        size_t perm_sz,
                        bool is_baldwinian);

int main(int argc, char **argv) {
        const double bin_cap = 20;
        const double prob_inst[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        const size_t perm1[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
        const size_t perm2[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        const size_t perm_sz = sizeof(perm1)/sizeof(*perm1);
        struct chromosome chrom1;
        chrom_init(&chrom1, false);
        struct bald_chrom chrom2;
        chrom_init(&chrom2.chrom, true);

        chrom1.perm = malloc(perm_sz * sizeof(*chrom1.perm));
        chrom2.chrom.perm = malloc(perm_sz * sizeof(*chrom2.chrom.perm));
        memcpy(chrom1.perm, perm1, perm_sz * sizeof(*chrom1.perm));
        memcpy(chrom2.chrom.perm, perm2, perm_sz * sizeof(*chrom1.perm));
        chrom_eval(&chrom1, prob_inst, perm_sz, bin_cap);
        chrom_eval(&chrom2.chrom, prob_inst, perm_sz, bin_cap);
        solution_copy(&chrom2.bald_sol, chrom2.chrom.sol);

        printf("chrom1:\n");
        chrom_print(&chrom1, perm_sz, false);
        printf("\nchrom2:\n");
        chrom_print(&chrom2.chrom, perm_sz, true);
        putchar('\n');

        for (int i = 0; i < 5; i++) {
                printf("mutating chrom1\n");
                chrom_mut(&chrom1, perm_sz);
                printf("evaluating chrom1\n");
                chrom_eval(&chrom1, prob_inst, perm_sz, bin_cap);
                printf("chrom1:\n");
                chrom_print(&chrom1, perm_sz, false);
        }
        putchar('\n');

        printf("crossover\n");
        struct chromosome child = chrom_cx(chrom1, chrom2.chrom, perm_sz);
        printf("evaluating child\n");
        chrom_eval(&child, prob_inst, perm_sz, bin_cap);
        printf("child:\n");
        chrom_print(&child, perm_sz, false);
        putchar('\n');

        printf("greedy lamarckian swap local search of child\n");
        printf("number of searches conducted: %d\n",
               chrom_search(&child, false, prob_inst, perm_sz, bin_cap,
                            true, SEARCHES, chrom_search_swap));
        printf("child:\n");
        chrom_print(&child, perm_sz, false);
        putchar('\n');

        printf("chrom2:\n");
        chrom_print(&chrom2.chrom, perm_sz, true);
        printf("greedy baldwinian swap local search of chrom2\n");
        printf("number of searches conducted: %d\n",
               chrom_search(&chrom2.chrom, true, prob_inst, perm_sz, bin_cap,
                            true, SEARCHES, chrom_search_swap));
        printf("chrom2:\n");
        chrom_print(&chrom2.chrom, perm_sz, true);
        putchar('\n');

        printf("chrom1:\n");
        chrom_print(&chrom1, perm_sz, false);
        printf("steep lamarckian swap local search of chrom1\n");
        printf("number of searches conducted: %d\n",
               chrom_search(&chrom1, false, prob_inst, perm_sz, bin_cap,
                            false, SEARCHES, chrom_search_swap));
        printf("chrom1:\n");
        chrom_print(&chrom1, perm_sz, false);
        putchar('\n');

        printf("chrom2:\n");
        chrom_print(&chrom2.chrom, perm_sz, true);
        printf("greedy baldwinian shuffle local search of chrom2\n");
        printf("number of searches conducted: %d\n",
               chrom_search(&chrom2.chrom, true, prob_inst, perm_sz, bin_cap,
                            true, SEARCHES, chrom_search_shuffle));
        printf("chrom2:\n");
        chrom_print(&chrom2.chrom, perm_sz, true);
        putchar('\n');

        chrom_destroy(&chrom1, false);
        chrom_destroy(&chrom2.chrom, true);
        chrom_destroy(&child, false);
        return 0;
}

static void chrom_print(const struct chromosome *chrom,
                        size_t perm_sz,
                        bool is_baldwinian) {
        printf("fitness: %lf\n"
               "permutation:\n",
               chrom->fitness);
        for (size_t i = 0; i < perm_sz; i++) {
                printf("%zu ", chrom->perm[i]);
        }
        putchar('\n');
        printf("solution:\n");
        solution_print(chrom->sol, stdout);
        if (is_baldwinian) {
                printf("baldwinian solution:\n");
                solution_print(((struct bald_chrom *)chrom)->bald_sol, stdout);
        }
}
