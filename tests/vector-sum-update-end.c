
/* Check that the vector_sum_update_start function works as advertised. */

#include <stdlib.h>
#include <stdio.h>

#include "../src/peakolator.c"

/* Probability a value is zero. This should be high so some actual compression
 * occurs. */
const double pr_zero = 0.99;

/* Maximum value. */
const val_t max_val = 1000;

/* Size of vector to generate. */
const size_t n = 500000;

/* Number of random updates to test. */
const size_t m = 10000;


/* Generate low-quality random values. */
val_t randval()
{
    if ((double) rand() / (double) RAND_MAX < pr_zero) {
        return 0;
    }
    else {
        return (val_t) (rand() % (int) max_val);
    }
}


void check_vector_sum(vector_sum_t* vs)
{
    if (vs->sum != vector_sum(vs->vec, vs->start, vs->end)) {
        fprintf(stderr, "Incorrect vector_sum_t.\n");
        exit(EXIT_FAILURE);
    }
}


int main()
{
    srand(1234);

    val_t* xs = malloc(n * sizeof(val_t));
    if (xs == NULL) {
        fprintf(stderr, "Cannot allocate %zu bytes", n * sizeof(val_t));
        return 1;
    }

    size_t i;
    for (i = 0; i < n; ++i) {
        xs[i] = randval();
    }

    vector_t* vec = vector_create(xs, n);
    if (vec == NULL) {
        fprintf(stderr, "Failed to create vector.");
        return 1;
    }

    vector_sum_t vs;
    vector_sum_set(&vs, vec, 0, n - 1);
    check_vector_sum(&vs);

    for (i = 0; i < m; ++i) {
        vector_sum_update_end(&vs, rand() % n);
        check_vector_sum(&vs);
    }

    return EXIT_SUCCESS;
}
