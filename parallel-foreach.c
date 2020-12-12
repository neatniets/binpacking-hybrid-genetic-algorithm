#include "parallel-foreach.h"
#include <pthread.h>
#include <assert.h>
#include <stdlib.h>

struct common_args {
        void *array;
        size_t count;
        size_t sz;
        void *context;
        para_foreach_func func;
        int num_threads;
        volatile int term_cond;
};
struct iterate_args {
        struct common_args *common;
        size_t index;
};
static void *pthread_iterate(void *args);

int parallel_foreach(int num_threads,
                     void *array,
                     size_t count,
                     size_t sz,
                     void *context,
                     para_foreach_func func) {
        int return_code = 0;
        struct common_args common;
        struct iterate_args *args = NULL;
        pthread_t *thrd_ids = NULL;

        if ((num_threads < 1)
            || (array == NULL)
            || (count < 1)
            || (sz < 1)
            || (func == NULL)) {
                return_code = ERR_BAD_ARGS;
                goto exit;
        }

        if (count < num_threads) {
                num_threads = count;
        }

        common = (struct common_args) {.array = array,
                                       .count = count,
                                       .sz = sz,
                                       .context = context,
                                       .func = func,
                                       .num_threads = num_threads,
                                       .term_cond = 0};
        args = malloc(num_threads * sizeof(*args));
        if (args == NULL) {
                return_code = ERR_MALLOC_FAIL;
                goto exit;
        }
        for (int i = 0; i < num_threads; i++) {
                args[i] = (struct iterate_args){.common = &common,
                                                .index = i};
        }

        thrd_ids = malloc((num_threads - 1) * sizeof(*thrd_ids));
        if (thrd_ids == NULL) {
                return_code = ERR_MALLOC_FAIL;
                goto exit;
        }

        for (int i = 0; i < num_threads - 1; i++) {
                int err = pthread_create(thrd_ids + i, NULL,
                                         pthread_iterate, args + i);
                if (err != 0) {
                        return_code = ERR_PTHREAD_FAIL;
                        goto exit;
                }
        }
        pthread_iterate(args + num_threads - 1);

        for (int i = 0; i < num_threads - 1; i++) {
                int err = pthread_join(thrd_ids[i], NULL);
                if (err != 0) {
                        return_code = ERR_PTHREAD_FAIL;
                        goto exit;
                }
        }

        return_code = common.term_cond;

exit:
        if (args != NULL) {
                free(args);
        }
        if (args != NULL) {
                free(thrd_ids);
        }
        return return_code;
}

static void *pthread_iterate(void *args) {
        struct iterate_args *argsv = args;
        for (size_t i = argsv->index;
             i < argsv->common->count;
             i += argsv->common->num_threads) {
                if (argsv->common->term_cond < 0) {
                        return NULL;
                }
                int err = argsv->common->func(argsv->common->array
                                              + (i * argsv->common->sz),
                                              argsv->common->context);
                if (err < 0) {
                        argsv->common->term_cond = err;
                        return NULL;
                }
        }
        return NULL;
}
