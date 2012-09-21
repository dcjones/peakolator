
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "peakolator.h"


/* Some helpful things to have. */


/* Malloc and exit on failure. */
void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


/* Realloc and exit on failure. */
void* realloc_or_die(void* ptr, size_t n)
{
    void* p = realloc(ptr, n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


/* Sparse vectors */


/* In principle, Peakolator operates on a array of uint32_t values. For
 * efficiency, we devide this array into chunks of size BLOCK_SIZE  and
 * precompute sums, and zero compress (i.e., leave out blocks with values
 * summing to 0) */

#define BLOCK_SIZE 1024



static idx_t idxmin(idx_t a, idx_t b)
{
    return a < b ? a : b;
}


static idx_t idxmax(idx_t a, idx_t b)
{
    return a > b ? a : b;
}


typedef struct block_t_
{
    /* Where this block occurs in the sequence. */
    idx_t idx;

    /* Sum over xs. */
    val_t sum;

    /* Values for idx, ..., idx + BLOCK_SIZE - 1 in the genomic sequence */
    val_t xs[BLOCK_SIZE];
} block_t;


/* Sum a part of a block.
 *
 * Args:
 *   block: A block.
 *   i: Start (within-block) index.
 *   j: End (within-block) index.
 *
 * Return:
 *   The sum xs[i] + ... + xs[j].
 *
 */
val_t block_sum(const block_t* block, idx_t i, idx_t j)
{
    val_t sum = 0;
    while (i <= j) sum += block->xs[i++];
    return sum;
}


/* A sparse vector representation. */
struct vector_t_
{
    /* Length of the sequence. */
    size_t n;

    /* Number of blocks. */
    size_t m;

    /* Constitutive blocks. */
    block_t* blocks;
};


/* Operations on vectors. */


/* Create a new sparse vector from a dense vector.
 *
 * Args:
 *   data: A dense vector of values.
 *   n: Length of the vector.
 *
 * Returns:
 *   An new sparse vector containing the same values as data.
 *
 */
vector_t* vector_create(const val_t* data, size_t n)
{
    vector_t* vec = malloc_or_die(sizeof(vector_t));

    /* Figure out how many blocks we need. */
    size_t i; /* genomic index */
    size_t m; /* block count */
    val_t block_sum = 0;

    for (i = 0, m = 0; i < n; ++i) {
        if (i > 0 && i % BLOCK_SIZE == 0 && block_sum > 0) {
            ++m;
            block_sum = 0;
        }

        block_sum += data[i];
    }

    if (block_sum > 0) ++m;
    ++m;

    /* Allocate */
    vec->n = n;
    vec->m = m;
    vec->blocks = malloc_or_die(m * sizeof(block_t));
    memset(vec->blocks, 0, m * sizeof(block_t));

    /* Initialize */
    size_t j; /* block index */
    size_t k; /* within block index */
    vec->blocks[j].idx = 0;
    for (i = 0, j = 0, k = 0, block_sum = 0; i < n; ++i, ++k) {
        if (i > 0 && i % BLOCK_SIZE == 0) {
            if (block_sum > 0) {
                vec->blocks[j].sum = block_sum;
                ++j;
            }

            vec->blocks[j].idx = i;
            k = 0;
            block_sum = 0;
        }

        vec->blocks[j].xs[k] = data[i];
        block_sum += data[i];
    }

    if (block_sum > 0) {
        vec->blocks[j].sum = block_sum;
    }

    return vec;
}


/* Free an allocated sparse vector.
 *
 * Args:
 *   vec: A sparse vector to be freed.
 */
void vector_free(vector_t* vec)
{
    free(vec->blocks);
    free(vec);
}


/* Find the block index corresponding to the given genomic index.
 *
 * Specifically, find the largest j, such that vec->blocks[j].idx <= i,
 * or, if there is no such j, the smallest j such that * vec->blocks[j].idx > i.
 *
 * Args:
 *   vec: A sparse vector.
 *   i: A genomic position.
 *   u: Lower bound for j.
 *   v: Upper bound for j;
 *
 * Returns:
 *   An block index j corresponding to the genomic index.
 *
 *
 */
idx_t vector_find_block_bound(const vector_t* vec, idx_t i, idx_t u, idx_t v)
{
    idx_t mid;
    while (u + 1 < v) {
        mid = u + (v - u) / 2;
        if (vec->blocks[mid].idx <= i) u = mid;
        else                           v = mid;
    }

    return u;
}


/* Find the largest block index j, such that vec->blocks[j].idx <= i.
 *
 * Args:
 *   vec: A sparse vector.
 *   i: A genomic position.
 *
 * Returns:
 *   An block index j corresponding to the genomic index.
 */
idx_t vector_find_block(const vector_t* vec, idx_t i)
{
    vector_find_block_bound(vec, i, 0, vec->m);
}


/* Find the sume of the values in the genomic interval [i, j].
 *
 * Args:
 *   vec: A sparse vector.
 *   i: Interval start.
 *   j: Interval end.
 *   u: Lower bound on block interval containing i.
 *   v: Upper bound on block interval containing i.
 *
 * Return:
 *   A sum of the values in [i, j].
 *
 */
val_t vector_sum_bound(const vector_t* vec, idx_t i, idx_t j, idx_t u, idx_t v)
{
    idx_t w = vector_find_block_bound(vec, i, u, v);

    val_t sum = 0;

    if (vec->blocks[w].idx <= j && vec->blocks[w].idx + BLOCK_SIZE - 1 >= i) {
        idx_t a = vec->blocks[w].idx >= i ? 0 : i - vec->blocks[w].idx;
        idx_t b = idxmin(j - vec->blocks[w].idx, BLOCK_SIZE - 1);
        sum += block_sum(&vec->blocks[w], a, b);
    }
    ++w;

    while (w < vec->m && vec->blocks[w].idx + BLOCK_SIZE - 1 <= j) {
        sum += vec->blocks[w].sum;
        ++w;
    }

    if (w < vec->m && vec->blocks[w].idx <= j) {
        sum += block_sum(&vec->blocks[w], 0, j - vec->blocks[w].idx);
    }

    return sum;
}


val_t vector_sum(const vector_t* vec, idx_t i, idx_t j)
{
    return vector_sum_bound(vec, i, j, 0, vec->m);
}

