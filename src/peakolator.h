
#ifndef PEAKOLATOR_H
#define PEAKOLATOR_H

#include <stdlib.h>
#include <stdint.h>

/* Zero compressed uint zectors. Peakolator scans across vectors cooresponding
 * to genomic sequences, with each position assigned a natural number.
 * */
typedef struct vector_t_ vector_t;


/* Vectors a sequence of values of type val_t, indexed with idx_t. */
typedef uint32_t val_t;
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
val_t vector_sum_bound(const vector_t* vec, idx_t i, idx_t j, idx_t u, idx_t v);


#endif

