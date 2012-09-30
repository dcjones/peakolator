
#ifndef PEAKOLATOR_H
#define PEAKOLATOR_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

/* Zero compressed vectors. Peakolator scans across vectors cooresponding
 * to genomic sequences, with each position assigned a natural number.
 * */
typedef struct vector_t_ vector_t;


/* Vectors a sequence of values of type val_t, indexed with idx_t. */
typedef double   val_t;
typedef uint32_t idx_t;


/* Create a new sparse vector from a dense vector.
 *
 * Args:
 *   data: A dense vector of values.
 *   n: Length of the vector.
 *
 * Returns:
 *   An new sparse vector containing the same values as data.
 */
vector_t* vector_create(const val_t* data, size_t n);


/* Free an allocated sparse vector.
 *
 * Args:
 *   vec: A sparse vector to be freed.
 */
void vector_free(vector_t* vec);


/* Length of the vector. */
idx_t vector_len(const vector_t* vec);


/* Find the sum of the values in genomic interval [i, j].
 *
 * Args:
 *   vec: A sparse vector.
 *   i: Interval start.
 *   j: Interval end.
 *
 * Returns:
 *   A sum of the values in [i, j].
 */
val_t vector_sum(const vector_t* vec, idx_t i, idx_t j);


/* Find the sum of the values in the genomic interval [i, j].
 *
 * This is more efficent than vector_sum, when it can be assumed that the intial
 * block index lies in [u, v].
 *
 * Args:
 *   vec: A sparse vector.
 *   i: Interval start.
 *   j: Interval end.
 *   u: Lower bound on block interval containing i.
 *   v: Upper bound on block interval containing i.
 *
 * Returns:
 *   A sum of the values in [i, j].
 */
val_t vector_sum_bound(const vector_t* vec, idx_t i, idx_t j, idx_t u, idx_t v,
                       idx_t* w0, idx_t* w1);


/* A representation of the interval [start, end] with an associated density. */
typedef struct interval_t_
{
    idx_t start, end;
    double density;
} interval_t;


/* Sort an array of interval in ascending order of density.
 *
 * Args:
 *   xs: An array of intervals.
 *   n: Number of intervals in xs.
 * */
void sort_intervals_asc_density(interval_t* xs, size_t n);


/* Sort an array of interval in descending order of density.
 *
 * Args:
 *   xs: An array of intervals.
 *   n: Number of intervals in xs.
 * */
void sort_intervals_des_density(interval_t* xs, size_t n);


/* Compare two intervals.
 *
 * Args:
 *   a: An interval.
 *   b: Another interval.
 *
 * Returns:
 *   Returns 0 if the intervals are equal (including equal densities),
 *   and non-zero otherwise. (Hint: this is just a wrapper around memcmp.)
 */
int interval_cmp(const interval_t* a, const interval_t* b);


/* A density function. */
typedef double (*density_function_t)(val_t, idx_t);


/* A prior over interval length. */
typedef double (*prior_function_t)(idx_t);


/* Run the peakolator algorithm, performing greedly one-dimensional clustering.
 *
 * The algorithm will repeatedly find the interval [i, j] in vec that maximized
 * f(x, k) + g(k), where x = vec[i] + ... + vec[j] and k = j - i + 1. This
 * higest-density interval is then subtracted from the search space and a new
 * highest-density interval is found. This process is repeated until no
 * intervals exist with f(x, k) + g(k) > -INFINITY.
 *
 * Typically (but not necessarily), f and g return log-probalilities.
 *
 * The function f must have the following properties:
 *   f(x', k) <= f(x, k) for all x' <= x
 *   f(x, k') <= f(x, k) for all k' >= k
 *
 * If there properties are not met, the algorithm will still run, but cannot
 * guaranteed to return an actual greedy clustering. (The highest-density
 * interval found at each step may be suboptimal.)
 *
 * Args:
 *   vec: Data to segment.
 *   f: Density function. (Behavior is undefined if this function is not a
 *      proper density function.)
 *   g: Prior on the length of intervals, or NULL for a flat prior.
 *   min_len: Minimum length of high-density intervals.
 *   max_len: Maximum length of high-density intervals.
 *   min_density: The minimum density of any high-density region reported.
 *   num_threads: Number of threads to use. If 0, use one thread per cpu core.
 *   out: A pointer to a pointer which will be set to an array holding the
 *        results of the clustering: high-density intervals in descending order
 *        of density.
 *
 * Returns:
 *   The number of high density intervals found.
 */
size_t peakolate(const vector_t* vec,
                 density_function_t f,
                 prior_function_t g,
                 idx_t min_len,
                 idx_t max_len,
                 double min_density,
                 unsigned int num_threads,
                 interval_t** out);

#endif

