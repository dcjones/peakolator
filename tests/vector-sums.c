
/* Check that sums computed on sparse vectors are correct. A large random
 * vector is generated "sparsified", and random sums are computed and checked.
 */

#include <stdlib.h>
#include <stdio.h>

#include "../src/peakolator.c"

/* Probability a value is zero. This should be high so some actual compression
 * occurs. */
const double pr_zero = 0.99;

/* Maximum value. */
const val_t max_val = 1000;

/* Size of vector to generate. */
const size_t n = 50000;

/* Number of random intervals to test. */
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


void randinterval(idx_t n, idx_t* u, idx_t* v)
{
    *u = rand() % n;
    *v = *u + (rand() % (n - *u));
}


val_t dense_sum(val_t* xs, idx_t i, idx_t j)
{
    val_t sum = 0;
    while (i <= j) sum += xs[i++];
    return sum;
}


#if 0
void dense_print(val_t* xs, idx_t i, idx_t j)
{
    while (i <= j) fprintf(stderr, " %u", xs[i++]);
    fprintf(stderr, "\n");
}
#endif


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

    idx_t u, v;
    for (i = 0; i < m; ++i) {
        randinterval(n, &u, &v);
        val_t true_sum   = dense_sum(xs, u, v);
        val_t sparse_sum = vector_sum(vec, u, v);

        if (true_sum != sparse_sum) {
            fprintf(stderr, "Incorrect sum.\n");
            return 1;
        }
    }

    vector_free(vec);
    free(xs);

    return 0;
}
