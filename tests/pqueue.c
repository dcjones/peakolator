
/* Test a priority queue with only one heap. */

#include <stdio.h>
#include <stdlib.h>

#include "../src/peakolator.c"


/* Comparison function for sorting interval bounds by ascending density_max */
static int interval_bound_cmp_asc(const void* a_, const void* b_)
{
    const interval_bound_t* a = (interval_bound_t*) a_;
    const interval_bound_t* b = (interval_bound_t*) b_;

    return -interval_bound_cmp_priority(a, b);
}


/* Sort an array of interval bounds in ascending order by density_max. */
static void sort_interval_bounds_asc(interval_bound_t* xs, size_t n)
{
    qsort(xs, n, sizeof(interval_bound_t), interval_bound_cmp_asc);
}


int main()
{
    pqueue_t* q = pqueue_create();
    if (q == NULL) {
        fprintf(stderr, "Failed to create a pqueue.\n");
        return EXIT_FAILURE;
    }

    size_t n = 100000;
    interval_bound_t* xs = malloc(n * sizeof(interval_bound_t));
    memset(xs, 0, n * sizeof(interval_bound_t));
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[i].start_min = 100000.0 * ((double) rand() / (double) RAND_MAX);
        xs[i].end_max = xs[i].start_min + 1000.0 * ((double) rand() / (double) RAND_MAX);
        xs[i].x_max = (double) rand() / (double) RAND_MAX;
        pqueue_enqueue(q, &xs[i]);
    }

    sort_interval_bounds_asc(xs, n);

    interval_bound_t x;
    for (i = 0; i < n; ++i) {
        if (!pqueue_dequeue(q, &x)) {
            fprintf(stderr, "Pqueue was prematurely empty.\n");
            return EXIT_FAILURE;
        }

        if (memcmp(&xs[i], &x, sizeof(interval_bound_t)) != 0) {
            fprintf(stderr, "Pop number %zu was incorrect.\n", i);
            return EXIT_FAILURE;
        }
    }

    pqueue_free(q);

    return EXIT_SUCCESS;
}

