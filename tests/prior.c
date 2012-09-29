
/* Test that priors over high-density interval length work correctly by
 * specifying an extremely strict prior that only allows intervals of a fixed
 * length. */

#include "../src/peakolator.c"

double f(val_t x, idx_t k)
{
    return log((double) x / (double) k);
}


const idx_t fixed_length = 50;

double g(idx_t k)
{
    return k == fixed_length ? 0.0 : -INFINITY;
}


int main()
{
    srand(1234);
    const size_t n = 1000;
    val_t* xs = malloc_or_die(n * sizeof(val_t));

    /* haystack */
    size_t i;
    for (i = 0; i < n; ++i) {
        xs[i] = rand() % 1000;
    }

    vector_t* vec = vector_create(xs, n);
    free(xs);

    interval_t* out;
    size_t count = peakolate(vec, f, g, 1, 10000, -INFINITY, 0, &out);

    if (count == 0) {
        fprintf(stderr, "No high density intervals were found.\n");
        return EXIT_FAILURE;
    }

    for (i = 0; i < count; ++i) {
        if (out[i].end - out[i].start + 1 != fixed_length) {
            fprintf(stderr,
                    "High-density interval of incorrect length reported.\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

