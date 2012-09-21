
#ifndef PEAKOLATOR_H
#define PEAKOLATOR_H

#include <stdlib.h>
#include <stdint.h>

/* TODO:
 *   We need a zero-compressed uint vector that supports reasonably efficient
 *   random access.
 *
 *   Similarly, a zero compressed binary vector.
 *   On second thought, let's just support the generalization that there is a
 *   maximum value at each position.
 *
 *   The interface will then consist of a function called "peakolate" that takes
 *   the compressed vector (or set of compressed vectors) and a callback score
 *   function (or a score function and a prior on length), as well as a callback
 *   function to enque new regions.
 */



/* Zero compressed uint zectors. Peakolator scans across vectors cooresponding
 * to genomic sequences, with each position assigned a natural number.
 * */
typedef struct vector_t_ vector_t;

typedef uint32_t val_t;
typedef uint32_t idx_t;

vector_t* vector_create(const val_t* data, size_t n);
void vector_free(vector_t* vec);
idx_t vector_find_block(const vector_t* vec, idx_t i);
val_t vector_sum_bound(const vector_t* vec, idx_t i, idx_t j, idx_t u, idx_t v);
val_t vector_sum(const vector_t* vec, idx_t i, idx_t j);

#endif

