
/* A needle-in-haystack test (find a single peak) in which the haystack is
 * all-zeros. */

#include <stdlib.h>
#include <math.h>

#include "../src/peakolator.c"


double f(val_t x, idx_t k)
{
    return log((double) x / (double) k);
}


double g(idx_t k __attribute__((unused)))
{
    return 0.0;
}


int main()
{
    srand(1234);
    const size_t n = 100000;
    val_t* xs = malloc_or_die(n * sizeof(val_t));

    /* haystack */
    memset(xs, 0, n * sizeof(val_t));

    /* needle */
    size_t i = rand() % n;
    xs[i] = 1000;

    vector_t* vec = vector_create(xs, n);
    free(xs);

    interval_t* out;
    size_t count = peakolate(vec, f, g, 1, 10000, -INFINITY, 0, &out);

    if (count == 0) {
        fprintf(stderr, "No high density intervals were found.\n");
        return EXIT_FAILURE;
    }
    else if (count > 1) {
        fprintf(stderr, "More than one high density interval was found.\n");
        return EXIT_FAILURE;
    }
    else if (out[0].start != i || out[0].end != i) {
        fprintf(stderr, "Incorrect needle found.\n");
        return EXIT_FAILURE;
    }

    vector_free(vec);
    free(out);

    return EXIT_SUCCESS;
}

