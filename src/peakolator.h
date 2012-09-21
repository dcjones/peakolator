
#ifndef PEAKOLATOR_H
#define PEAKOLATOR_H

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





#endif

