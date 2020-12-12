#ifndef PARALLEL_FOREACH_H
#define PARALLEL_FOREACH_H

#include <stddef.h>

/** Will return negative non-zero if all threads should terminate;
  * returning >=0 will continue thread */
typedef int (*para_foreach_func)(void *elem,
                                 void *context);

enum parallel_foreach_errors {
        ERR_BAD_ARGS = 1,
        ERR_MALLOC_FAIL,
        ERR_PTHREAD_FAIL
};
/** Returns 0 if all threads ran to completion, else returns the negative
 * non-zero value that terminated the threads or an error code for internal
 * errors */
int parallel_foreach(int num_threads,
                     void *array,
                     size_t count,
                     size_t sz,
                     void *context,
                     para_foreach_func func);

#endif /* !PARALLEL_FOREACH_H */
