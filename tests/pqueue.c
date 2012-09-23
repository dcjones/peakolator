
#include <stdio.h>
#include <stdlib.h>

#include "../src/peakolator.h"

int main()
{
    pqueue_t* q = pqueue_create();
    if (q == NULL) {
        fprintf(stderr, "Failed to create a pqueue.\n");
        return EXIT_FAILURE;
    }

    size_t n = 100000;
    interval_t* xs = malloc(n * sizeof(interval_t));
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[i].density = (double) rand() / (double) RAND_MAX;
        pqueue_enqueue(q, &xs[i]);
    }

    sort_intervals_asc(xs, n);

    interval_t x;
    for (i = 0; i < n; ++i) {
        if (!pqueue_dequeue(q, &x)) {
            fprintf(stderr, "Pqueue was prematurely empty.\n");
            return EXIT_FAILURE;
        }

        if (interval_cmp(&xs[i], &x) != 0) {
            fprintf(stderr, "Pop number %zu was incorrect.\n", i);
            return EXIT_FAILURE;
        }
    }

    pqueue_free(q);

    return EXIT_SUCCESS;
}

