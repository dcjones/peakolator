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


void pthread_mutex_init_or_die(pthread_mutex_t* mutex,
                               const pthread_mutexattr_t* attr)
{
    if (pthread_mutex_init(mutex, attr) != 0) {
        fprintf(stderr, "Failed to init mutex.\n");
        exit(EXIT_FAILURE);
    }
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


#if 0
static idx_t idxmax(idx_t a, idx_t b)
{
    return a > b ? a : b;
}
#endif


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


/* Create a new sparse vector from a dense vector. */
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


/* Free an allocated sparse vector. */
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
    return vector_find_block_bound(vec, i, 0, vec->m);
}


/* Find the sum of the values in the genomic interval [i, j]. */
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


/* Find the sum of the values in the genomic interval [i, j]. */
val_t vector_sum(const vector_t* vec, idx_t i, idx_t j)
{
    return vector_sum_bound(vec, i, j, 0, vec->m);
}


/* Copy an interval. */
static void interval_copy(interval_t* dest, const interval_t* src)
{
    memcpy(dest, src, sizeof(interval_t));
}


/* Swap two intervals. */
static void interval_swap(interval_t* a, interval_t* b)
{
    interval_t tmp;
    interval_copy(&tmp, a);
    interval_copy(a, b);
    interval_copy(b, &tmp);
}


/* Comparison function for sorting intervals in asceding order. */
static int interval_cmp_asc(const void* a_, const void* b_)
{
    const interval_t* a = (interval_t*) a_;
    const interval_t* b = (interval_t*) b_;

    if      (a->density == b->density) return 0;
    else if (a->density <  b->density) return 1;
    else                               return -1;
}


/* Sort an array of interval in descending order. */
void sort_intervals_asc(interval_t* xs, size_t n)
{
    qsort(xs, n, sizeof(interval_t), interval_cmp_asc);
}


/* Compare two intervals. */
int interval_cmp(const interval_t* a, const interval_t* b)
{
    return memcmp(a, b, sizeof(interval_t));
}


/* A concurrent priority queue for genomic intervals.
 *
 * This is just a coarsely locked binary max-heap, nothing fancy.
 */
struct pqueue_t_
{
    /* Coarse lock. */
    pthread_mutex_t mutex;

    /* Number of items in the heap. */
    size_t n;

    /* Size reserved. */
    size_t size;

    /* Items. */
    interval_t* xs;
};


/* Create a new priority queue. */
pqueue_t* pqueue_create()
{
    pqueue_t* q = malloc_or_die(sizeof(pqueue_t));
    q->n = 0;
    q->size = 1024;
    q->xs = malloc_or_die(q->size * sizeof(interval_t));
    pthread_mutex_init_or_die(&q->mutex, NULL);
    return q;
}


/* Free an priority queue created with pqueue_create. */
void pqueue_free(pqueue_t* q)
{
    pthread_mutex_destroy(&q->mutex);
    free(q->xs);
    free(q);
}


/* Double the size reserved for a heap. */
static void pqueue_expand(pqueue_t* q)
{
    q->size *= 2;
    q->xs = realloc_or_die(q->xs, q->size * sizeof(interval_t));
}


/* Binary heap index function. */
static inline size_t pqueue_parent_idx (size_t i) { return (i - 1) / 2; }
static inline size_t pqueue_left_idx   (size_t i) { return 2 * i + 1; }
static inline size_t pqueue_right_idx  (size_t i) { return 2 * i + 2; }


/* Enqueue an item. */
void pqueue_enqueue(pqueue_t* q, const interval_t* interval)
{
    pthread_mutex_lock(&q->mutex);

    /* insert */
    if (q->n == q->size) pqueue_expand(q);
    size_t j, i = q->n++;
    interval_copy(&q->xs[i], interval);

    /* percolate up */
    while (i > 0) {
        j = pqueue_parent_idx(i);
        if (q->xs[j].density < q->xs[i].density) {
            interval_swap(&q->xs[i], &q->xs[j]);
            i = j;
        }
        else break;
    }

    pthread_mutex_unlock(&q->mutex);
}


/* Pop the item with the largest density from the queue. */
bool pqueue_dequeue(pqueue_t* q, interval_t* interval)
{
    pthread_mutex_lock(&q->mutex);

    if (q->n == 0) {
        pthread_mutex_unlock(&q->mutex);
        return false;
    }

    interval_copy(interval, &q->xs[0]);

    /* replace head */
    interval_copy(&q->xs[0], &q->xs[--q->n]);

    /* percolate down */
    size_t l, r, j, i = 0;
    while (true) {
        l = pqueue_left_idx(i);
        r = pqueue_right_idx(i);

        if (l < q->n) {
            j = r < q->n && q->xs[r].density > q->xs[l].density ? r : l;

            if (q->xs[j].density > q->xs[i].density) {
                interval_swap(&q->xs[j], &q->xs[i]);
                i = j;
            }
            else break;
        }
        else break;
    }

    pthread_mutex_unlock(&q->mutex);
    return true;
}



/* Just sketching stuff out here. */
#if 0

static void* peakolator_thread(void* param)
{

}


void peakolate_async(const vector_t* vec,
                     density_function_t f,
                     prior_function_t g,
                     pqueue_t* out,
                     pthread_cond_t* cond)
{
    /* How many threads to launch? */
    size_t num_threads = 8;
    pthread_t* threads = malloc_or_die(num_threads * sizeof(pthread_t));

    size_t i;
    for (i = 0; i < num_threads; ++i) {
        /*pthread_create(&theads[i], NULL, peakolator_thread, XXX);*/
    }


    free(threads);
}


void peoklate(const vector_t* vec,
              density_function_t f,
              prior_function_t g)
{
    pqueue_t* out;
    /* TODO: create */

    /* TODO: Check the output of pthread functions. */
    pthread_cond_t cond;
    pthread_cond_init(&cond, NULL);

    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_lock(&mutex);

    peakolate_async(vec, f, g, out, &cond);

    pthread_cond_wait(&cond, &mutex);

    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&mutex);

    /* TODO: Return something. Dump `out`. */
}

#endif


