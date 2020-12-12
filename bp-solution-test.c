#include "bp-solution.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
        struct solution sol;
        solution_init(&sol);
        const double arr[] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
        const size_t perm[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        const size_t arr_sz = sizeof(arr)/sizeof(*arr);
        const double bin_cap = 20;

        solution_first_fit(&sol, arr, arr_sz, perm, bin_cap);
        solution_print(sol, stdout);
        printf("solution Falkenauer fitness: %lf\n",
               solution_eval(sol, bin_cap));

        size_t *new_perm = solution_reverse_first_fit(sol, arr_sz);
        for (size_t i = 0; i < arr_sz; i++) {
                printf("%zu ", new_perm[i]);
        }
        putchar('\n');

        solution_first_fit(&sol, arr, arr_sz, new_perm, bin_cap);
        solution_print(sol, stdout);
        printf("solution Falkenauer fitness: %lf\n",
               solution_eval(sol, bin_cap));

        free(new_perm);
        solution_destroy(sol);
        return 0;
}
