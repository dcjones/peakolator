#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "peakolator.h"


/* Some helpful things to have. */


/* Number of cpus. */
unsigned int ncpus()
{
    /* TODO: This is the only thing in the codebase that is not entirely
     * portable. We might want to at least give windows some love. */
    return sysconf(_SC_NPROCESSORS_ONLN);
}


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


/* Call pthread_mutex_init and exit on failure. */
void pthread_mutex_init_or_die(pthread_mutex_t* mutex,
                               const pthread_mutexattr_t* attr)
{
    if (pthread_mutex_init(mutex, attr) != 0) {
        fprintf(stderr, "Failed to init pthreads mutex.\n");
        exit(EXIT_FAILURE);
    }
}


/* Call pthread_cond_init and exit on failure. */
void pthread_cond_init_or_die(pthread_cond_t* cond,
                              const pthread_condattr_t* attr)
{
    if (pthread_cond_init(cond, attr) != 0) {
        fprintf(stderr, "Failed to init pthreads cond.\n");
        exit(EXIT_FAILURE);
    }
}


/* Sparse vectors */


/* In principle, Peakolator operates on a array of uint32_t values. For
 * efficiency, we devide this array into chunks of size BLOCK_SIZE  and
 * precompute sums, and zero compress (i.e., leave out blocks with values
 * summing to 0) */

#define BLOCK_SIZE 32


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
static val_t block_sum(const block_t* block, idx_t i, idx_t j)
{
    assert(i < BLOCK_SIZE);
    assert(j < BLOCK_SIZE);

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

    /* Allocate */
    vec->n = n;
    vec->m = m;
    vec->blocks = malloc_or_die(m * sizeof(block_t));
    memset(vec->blocks, 0, m * sizeof(block_t));

    /* Initialize */
    size_t j; /* block index */
    size_t k; /* within block index */
    for (i = 0, j = 0, k = 0, block_sum = 0; i < n; ++i, ++k) {
        if (i > 0 && i % BLOCK_SIZE == 0) {
            if (block_sum > 0) {
                vec->blocks[j].sum = block_sum;
                ++j;
                if (j >= m) break;
            }

            vec->blocks[j].idx = i;
            k = 0;
            block_sum = 0;
        }

        vec->blocks[j].xs[k] = data[i];
        block_sum += data[i];
    }

    if (j < m && block_sum > 0) {
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


/* Length of the vector. */
idx_t vector_len(const vector_t* vec)
{
    return vec->n;
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
 */
static idx_t vector_find_block_bound(const vector_t* vec, idx_t i, idx_t u, idx_t v)
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
val_t vector_sum_bound(const vector_t* vec, idx_t i, idx_t j, idx_t u, idx_t v,
                       idx_t* w0, idx_t* w1)
{
    idx_t w = vector_find_block_bound(vec, i, u, v);
    if (w0) *w0 = w;

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

    if (w1) {
        if (w > 0 && (w == vec->m || vec->blocks[w].idx > j)) --w;
        *w1 = w;
    }

    return sum;
}


/* Find the sum of the values in the genomic interval [i, j]. */
val_t vector_sum(const vector_t* vec, idx_t i, idx_t j)
{
    return vector_sum_bound(vec, i, j, 0, vec->m, NULL, NULL);
}


/* An interval with an associated sum that can be efficiently update.
 *
 * In particular, we are interested in updating the interval [i, j] to
 * [i + delta, j] or [i, j + delta] while doing the minimal work to recompute
 * the sum of that interval.
 */
typedef struct vector_sum_t_
{
    const vector_t* vec;
    idx_t start, end;
    idx_t start_block, end_block;
    val_t sum;
} vector_sum_t;


/* Create a num vector_sum_t. */
void vector_sum_set(vector_sum_t* vs, const vector_t* vec,
                    idx_t start, idx_t end)
{
    vs->vec = vec;
    vs->start = start;
    vs->end   = end;
    vs->sum   = vector_sum_bound(vec, start, end,
                                 0, vec->m,
                                 &vs->start_block, &vs->end_block);
}


/* Adjust the start of a vector_sum_t. */
void vector_sum_update_start(vector_sum_t* vs, idx_t new_start)
{
    assert(new_start <= vs->end);

    idx_t w0, w1;

    /* Subtracting from the sum. */
    if (new_start > vs->start) {
        val_t off = vector_sum_bound(vs->vec, vs->start, new_start - 1,
                                     vs->start_block, vs->start_block + 1,
                                     &w0, &w1);
        vs->start = new_start;
        vs->start_block = w1;
        vs->sum -= off;
    }
    /* Adding to the sum. */
    else if (new_start < vs->start) {
        val_t off = vector_sum_bound(vs->vec, new_start, vs->start - 1,
                                     0, vs->start_block + 1,
                                     &w0, &w1);
        vs->start = new_start;
        vs->start_block = w0;
        vs->sum += off;
    }

    if (vs->start_block + 1 < vs->vec->m &&
        vs->start >= vs->vec->blocks[vs->start_block + 1].idx) {
        ++vs->start_block;
    }
    else if (vs->start_block > 0 &&
             vs->start < vs->vec->blocks[vs->start_block].idx) {
        --vs->end_block;
    }
}


void vector_sum_inc_start(vector_sum_t* vs)
{
    assert(vs->start < vs->end);
    assert(vs->start_block < vs->vec->m);
    const block_t* block = &vs->vec->blocks[vs->start_block];

    if (vs->start >= block->idx &&
        vs->start < block->idx + BLOCK_SIZE) {
        vs->sum -= block->xs[vs->start - block->idx];
    }

    if (vs->start < block->idx + BLOCK_SIZE &&
        vs->start_block + 1 < vs->vec->m &&
        vs->start + 1 >= vs->vec->blocks[vs->start_block + 1].idx) {
        ++vs->start_block;
    }

    ++vs->start;
}


/* Adjust the end of a vector_sum_t */
void vector_sum_update_end(vector_sum_t* vs, idx_t new_end)
{
    assert(new_end >= vs->start);

    idx_t w0, w1;

    /* Adding to sum. */
    if (new_end > vs->end) {
        val_t off = vector_sum_bound(vs->vec, vs->end + 1, new_end,
                                     vs->end_block, vs->end_block + 1,
                                     &w0, &w1);
        vs->end = new_end;
        vs->end_block = w1;
        vs->sum += off;
    }
    /* Subtracting from sum. */
    else if (new_end < vs->end) {
        val_t off = vector_sum_bound(vs->vec, new_end + 1, vs->end,
                                     vs->start_block, vs->end_block + 1,
                                     &w0, &w1);
        vs->end = new_end;
        vs->end_block = w0;
        vs->sum -= off;
    }

    if (vs->end_block + 1 < vs->vec->m &&
        vs->end >= vs->vec->blocks[vs->end_block + 1].idx) {
        ++vs->end_block;
    }
    else if (vs->end_block > 0 &&
             vs->end < vs->vec->blocks[vs->end_block].idx) {
        --vs->end_block;
    }
}


void vector_sum_inc_end(vector_sum_t* vs)
{
    if (vs->end < vs->vec->blocks[vs->end_block].idx + BLOCK_SIZE &&
        vs->end_block +1 < vs->vec->m &&
        vs->end + 1 >= vs->vec->blocks[vs->end_block + 1].idx) {
        ++vs->end_block;
    }
    ++vs->end;

    const block_t* block = &vs->vec->blocks[vs->end_block];

    if (vs->end >= block->idx &&
        vs->end <  block->idx + BLOCK_SIZE) {
        vs->sum += block->xs[vs->end - block->idx];
    }
}


/* A bound on the start and end of an interval. */
typedef struct interval_bound_t_
{
    /* Inclusive interval in which in which the start in bound. */
    idx_t start_min, start_max;

    /* Inclusive interval in which in which the end in bound. */
    idx_t end_min, end_max;

    /* An upper bound on the density of any interval within this bound. */
    double density_max;

    /* An upper bound on the mass of any interval within this bound. */
    val_t x_max;
} interval_bound_t;


/* Minimum interval length within this bound. */
static idx_t interval_bound_min_len(const interval_bound_t* bound)
{
    if (bound->end_min > bound->start_max) {
        return bound->end_min - bound->start_max + 1;
    }
    else return 1;
}


/* Compare two interval bounds using heuristic ordering.
 *
 * Specifically, the bound with the largest upper-bound density will be first,
 * with naive density (mass/area) used as a tie breker.
 *
 * Args:
 *   a: An interval bound.
 *   b: An interval bound.
 *
 * Returns:
 *   > 0 if a is higher priority than b, 0 if equal, and < 0 if lesser.
 */
static int interval_bound_cmp_priority(const interval_bound_t* a,
                                       const interval_bound_t* b)
{
    if      (a->density_max > b->density_max) return  1;
    else if (a->density_max < b->density_max) return -1;
    else {
        double u = (double) a->x_max / (double) (a->end_max - a->start_min + 1);
        double v = (double) b->x_max / (double) (b->end_max - b->start_min + 1);

        if (u > v) return  1;
        if (u < v) return -1;
        else       return  0;
    }
}


/* Arithmetic with uint64_t that bottoms out rather than overflows. */
static uint64_t uint64_add(uint64_t a, uint64_t b)
{
    return a < UINT64_MAX - b ? a + b : UINT64_MAX;
}


static uint64_t uint64_sub(uint64_t a, uint64_t b)
{
    return a > b ? a - b : 0;
}


static uint64_t uint64_mul(uint64_t a, uint64_t b)
{
    return b == 0 || a <= UINT64_MAX / b ? a * b : UINT64_MAX;
}



/* Count the number of subintervals of length k contained in an interval of size
 * n, for all k in [k_min, k_max].
 *
 * Args:
 *   n: Size of containing interval.
 *   k_min: Minimum subinterval length.
 *   k_max: Maximum subinterval length.
 *
 * Returns:
 *   Number of subintervals, or UINT64_MAX if there are too many to count.
 */
static uint64_t subinterval_count(idx_t n_, idx_t k_min_, idx_t k_max_)
{
    uint64_t n = n_;
    uint64_t k_min = k_min_ < n_ ? k_min_ : n_;
    uint64_t k_max = k_max_ < n_ ? k_max_ : n_;

    /* Two times number of intervals with length >= k_min. */
    uint64_t u = uint64_mul(n - k_min - 1, n - k_min);

    /* Two times number of intervals with length > k_max. */
    uint64_t v = uint64_mul(n - k_max, n - k_max + 1);

    return (u - v) / 2;
}


/* Count the number of intervals contained within the bound. */
static uint64_t interval_bound_count(const interval_bound_t* bound,
                                     idx_t k_min, idx_t k_max)
{
    /* Interval bounds in peakolator are always either equal or disjoint. */
    if (bound->start_min == bound->end_min &&
       bound->start_max == bound->end_max) {
        return subinterval_count(bound->end_max - bound->start_min + 1,
                                 k_min, k_max);
    }
    else if (bound->start_max < bound->end_min) {
        /* total */
        uint64_t t = subinterval_count(bound->end_max - bound->start_min + 1,
                                       k_min, k_max);

        /* left */
        uint64_t l = subinterval_count((bound->end_min - 1) -
                                       bound->start_min + 1,
                                       k_min, k_max);

        /* right */
        uint64_t r = subinterval_count(bound->end_max -
                                       (bound->start_max + 1) + 1,
                                       k_min, k_max);
        /* center */
        uint64_t c = subinterval_count((bound->end_min - 1) -
                                       (bound->start_max + 1) + 1,
                                       k_min, k_max);

        return uint64_sub(uint64_sub(uint64_add(t, c), l), r);
    }
    else {
        /* It's not really zero, but it's tricky and we don't need to know. */
        return 0;
    }
}


/* Copy an interval bound. */
static void interval_bound_copy(interval_bound_t* dest,
                                const interval_bound_t* src)
{
    dest->start_min   = src->start_min;
    dest->start_max   = src->start_max;
    dest->end_min     = src->end_min;
    dest->end_max     = src->end_max;
    dest->density_max = src->density_max;
    dest->x_max       = src->x_max;
}


/* Swap two interval bounds. */
static void interval_bound_swap(interval_bound_t* a,
                                interval_bound_t* b)
{
    interval_bound_t tmp;
    interval_bound_copy(&tmp, a);
    interval_bound_copy(a, b);
    interval_bound_copy(b, &tmp);
}


/* Copy an interval. */
static void interval_copy(interval_t* dest, const interval_t* src)
{
    dest->start   = src->start;
    dest->end     = src->end;
    dest->density = src->density;
    dest->x       = src->x;
}


/* Swap two intervals. */
#if 0
static void interval_swap(interval_t* a, interval_t* b)
{
    interval_t tmp;
    interval_copy(&tmp, a);
    interval_copy(a, b);
    interval_copy(b, &tmp);
}
#endif


/* Comparison function for sorting intervals in asceding order. */
static int interval_cmp_asc_density(const void* a_, const void* b_)
{
    const interval_t* a = (interval_t*) a_;
    const interval_t* b = (interval_t*) b_;

    if      (a->density == b->density) return 0;
    else if (a->density <  b->density) return -1;
    else                               return 1;
}


/* Comparison function for sorting intervals in asceding order. */
static int interval_cmp_des_density(const void* a_, const void* b_)
{
    const interval_t* a = (interval_t*) a_;
    const interval_t* b = (interval_t*) b_;

    if      (a->density == b->density) return 0;
    else if (a->density <  b->density) return 1;
    else                               return -1;
}


/* Comparison function for sorting intervals in ascending position. */
static int interval_cmp_asc_pos(const void* a_, const void* b_)
{
    const interval_t* a = (interval_t*) a_;
    const interval_t* b = (interval_t*) b_;
    if (a->start < b->start) return -1;
    if (a->start > b->start) return  1;
    else {
        if      (a->end < b->end) return -1;
        else if (a->end > b->end) return  1;
        else return 0;
    }
}


/* Sort an array of intervals by density in ascending order. */
void sort_intervals_asc_density(interval_t* xs, size_t n)
{
    qsort(xs, n, sizeof(interval_t), interval_cmp_asc_density);
}


/* Sort an array of intervals by density in descending order. */
void sort_intervals_des_density(interval_t* xs, size_t n)
{
    qsort(xs, n, sizeof(interval_t), interval_cmp_des_density);
}


/* Sort intervals by position. */
void sort_intervals_asc_pos(interval_t* xs, size_t n)
{
    qsort(xs, n, sizeof(interval_t), interval_cmp_asc_pos);
}



/* Compare two intervals. */
int interval_cmp(const interval_t* a, const interval_t* b)
{
    return memcmp(a, b, sizeof(interval_t));
}


/* A stack of intervals. */
typedef struct interval_stack_t_
{
    /* Items in the stack. */
    interval_t* xs;

    /* Number of items. */
    size_t n;

    /* Number of pending items. */
    size_t pending;

    /* Size reserved. */
    size_t size;

    pthread_mutex_t mutex;
    pthread_cond_t cond;

} interval_stack_t;


/* Create an interval stack.
 *
 * Returns:
 *   An interval stack.
 */
static interval_stack_t* interval_stack_create()
{
    interval_stack_t* s = malloc_or_die(sizeof(interval_stack_t));
    s->n = 0;
    s->pending = 0;
    s->size = 1024;
    s->xs = malloc_or_die(s->size * sizeof(interval_t));
    pthread_mutex_init_or_die(&s->mutex, NULL);
    pthread_cond_init_or_die(&s->cond, NULL);
    return s;
}


/* Free a interval stack.
 *
 * Args:
 *   s: An interval stack.
 */
static void interval_stack_free(interval_stack_t* s)
{
    free(s->xs);
    pthread_mutex_destroy(&s->mutex);
    pthread_cond_destroy(&s->cond);
    free(s);
}


/* Double the size reserved for an interval stack. */
static void interval_stack_expand(interval_stack_t* s)
{
    s->size *= 2;
    s->xs = realloc_or_die(s->xs, s->size * sizeof(interval_t));
}


/* Push an interval onto an interval stack.
 *
 * Args:
 *   s: An interval stack.
 *   interval: An interval to push.
 */
static void interval_stack_push(interval_stack_t* s, const interval_t* interval)
{
    pthread_mutex_lock(&s->mutex);
    if (s->n == s->size) interval_stack_expand(s);
    interval_copy(&s->xs[s->n++], interval);
    pthread_mutex_unlock(&s->mutex);
}


/* Pop an item from an interval stack.
 *
 * Args:
 *   s: An interval stack.
 *   interval: A pointer where the popped item should be copied.
 *
 * Returns:
 *   true if an item was poped (false if the stack is empty).
 */
static bool interval_stack_pop(interval_stack_t* s, interval_t* interval)
{
    pthread_mutex_lock(&s->mutex);

    while (s->n == 0 && s->pending > 0) {
        pthread_cond_wait(&s->cond, &s->mutex);
    }

    if (s->n == 0) {
        pthread_mutex_unlock(&s->mutex);
        return false;
    }

    interval_copy(interval, &s->xs[--s->n]);
    ++s->pending;
    pthread_mutex_unlock(&s->mutex);
    pthread_cond_signal(&s->cond);
    return true;
}


/* Return true if the stack is finished.
 *
 * Finished in this context means that it is empty and there are no pending
 * items (i.e., interval_stack_finish_one was called on everything that was
 * dequeued).
 *
 * Args:
 *   s: A stack.
 *
 * Returns
 *   true if finished.
 */
static bool interval_stack_finished(interval_stack_t* s)
{
    pthread_mutex_lock(&s->mutex);
    bool finished = s->n == 0 && s->pending == 0;
    pthread_mutex_unlock(&s->mutex);
    return finished;
}


/* Notify the stack that a previously popped interval is no longer being
 * processed.
 *
 * This function needs to be called for every popped interval so the stack can
 * keep track of when no more intervals might be pushed.
 */
static void interval_stack_finish_one(interval_stack_t* s)
{
    pthread_mutex_lock(&s->mutex);
    --s->pending;
    pthread_mutex_unlock(&s->mutex);
    pthread_cond_broadcast(&s->cond);
}


/* A concurrent priority queue for genomic intervals.
 *
 * This is a partitioned priority queue. To reduce contention, there are k
 * independent queues each with their own lock.
 */
typedef struct pqueue_t_
{
    /* Number of items in the heap. */
    size_t n;

    /* Size reserved. */
    size_t size;

    /* Items. */
    interval_bound_t* xs;
} pqueue_t;


/* Create a new priority queue. */
static pqueue_t* pqueue_create()
{
    pqueue_t* q = malloc_or_die(sizeof(pqueue_t));
    q->n = 0;
    q->size = 1024;
    q->xs = malloc_or_die(q->size * sizeof(interval_bound_t));

    return q;
}


/* Free an priority queue created with pqueue_create. */
static void pqueue_free(pqueue_t* q)
{
    free(q->xs);
    free(q);
}


/* Double the size reserved for heap i. */
static void pqueue_expand(pqueue_t* q)
{
    q->size *= 2;
    q->xs = realloc_or_die(q->xs,
                           q->size * sizeof(interval_bound_t));
}


/* Binary heap index function. */
static inline size_t pqueue_parent_idx (size_t i) { return (i - 1) / 2; }
static inline size_t pqueue_left_idx   (size_t i) { return 2 * i + 1; }
static inline size_t pqueue_right_idx  (size_t i) { return 2 * i + 2; }


/* Insert an item into the priority queue in heap number i.
 *
 * Args:
 *   h: Heap number.
 *   q: A pqueue.
 *   interval: Interval to enqueue.
 */
static void pqueue_enqueue(pqueue_t* q, const interval_bound_t* bound)
{
    /* insert */
    if (q->n == q->size) pqueue_expand(q);

    size_t j, i = q->n++;
    interval_bound_copy(&q->xs[i], bound);

    /* percolate up */
    while (i > 0) {
        j = pqueue_parent_idx(i);
        if (interval_bound_cmp_priority(&q->xs[i], &q->xs[j]) > 0) {
            interval_bound_swap(&q->xs[i], &q->xs[j]);
            i = j;
        }
        else break;
    }
}


/* Pop the item with the largest density from the queue.
 *
 * If the queue is empty, this function will block until an item becomes
 * available, unless pqueue_finish has been called, in which case it will return
 * immediately.
 *
 * Args:
 *   h: Heap number.
 *   q: A pqueue.
 *   interval: Location to copy popped interval.
 *
 * Returns:
 *   false if the queue is finished and empty.
 */
static bool pqueue_dequeue(pqueue_t* q, interval_bound_t* bound)
{
    if (q->n == 0) return false;

    interval_bound_copy(bound, &q->xs[0]);
    if (--q->n == 0) return true;

    /* replace head */
    interval_bound_copy(&q->xs[0], &q->xs[q->n]);

    /* percolate down */
    size_t l, r, j, i = 0;
    while (true) {
        l = pqueue_left_idx(i);
        r = pqueue_right_idx(i);

        if (l >= q->n) break;

        if (r >= q->n) {
            j = l;
        }
        else {
            if (interval_bound_cmp_priority(&q->xs[l], &q->xs[r]) > 0) {
                j = l;
            }
            else {
                j = r;
            }
        }

        if (interval_bound_cmp_priority(&q->xs[i], &q->xs[j]) < 0) {
            interval_bound_swap(&q->xs[j], &q->xs[i]);
            i = j;
        }
        else break;
    }

    return true;
}


/* Clear a priority queue. */
static void pqueue_clear(pqueue_t* q)
{
    q->n = 0;
}


/* A lookup table for functions of the form idx_t -> double. */
typedef struct prior_lookup_t_
{
    idx_t min, max;
    double* xs;
} prior_lookup_t;


/* Create a lookup table from the given function.
 *
 * Args:
 *   g: A prior function.
 *   min: Minimum value the lookup table will be evaluated at.
 *   max: Maximum value the lookup table will be evaluated at.
 *
 * Returns:
 *   A lookup table.
 */
static prior_lookup_t* memoize_prior(prior_function_t g, idx_t min, idx_t max)
{
    prior_lookup_t* lookup = malloc_or_die(sizeof(prior_lookup_t));
    lookup->xs = malloc_or_die((max - min + 1) * sizeof(double));
    lookup->min = min;
    lookup->max = max;

    idx_t i;
    for (i = min; i <= max; ++i) {
        lookup->xs[i - min] = g(i);
    }

    return lookup;
}


/* Build a lookup table for a prior functions mode.
 *
 * Specifically, for a given number l, we need to quickly compute the function
 *     max g(k) for k <= l,
 *
 * Args:
 *   g: A prior function.
 *   min: Minimum value the lookup table will be evaluated at.
 *   max: Maximum value the lookup table will be evaluated at.
 *
 * Returns:
 *   A lookup table.
 */
static prior_lookup_t* memoize_prior_mode(prior_function_t g,
                                          idx_t min, idx_t max)
{
    prior_lookup_t* lookup = malloc_or_die(sizeof(prior_lookup_t));
    lookup->xs = malloc_or_die((max - min + 1) * sizeof(double));
    lookup->min = min;
    lookup->max = max;

    double v, mode = -INFINITY;
    idx_t i;
    for (i = min; i <= max; ++i) {
        v = g(i);
        if (v > mode) mode = v;
        lookup->xs[i - min] = mode;
    }

    return lookup;
}


/* Evaluate the function represented by a lookup table.
 *
 * Args:
 *   lookup: Lookup table.
 *   k: Value at which to evaluate the function.
 *
 * Returns:
 *   The value of the function at k, or if k < k_min or k > k_max, -INFINITY.
 */
static double prior_lookup_eval(const prior_lookup_t* lookup,
                                  idx_t k)
{
    if (k < lookup->min) return -INFINITY;
    else if (k > lookup->max) k = lookup->max;
    return lookup->xs[k - lookup->min];
}


/* Free a lookup table. */
static void prior_lookup_free(prior_lookup_t* lookup)
{
    free(lookup->xs);
    free(lookup);
}


/* Compute an upper bound on the maximum density interval.
 *
 * Args:
 *   f: Density function.
 *   k_min: The smallest interval considered.
 *   g_mode: A lookup table for prior function mode.
 *   x: Total mass in the interval.
 *   k: Length of the interval.
 *
 * Returns:
 *   An upper bound density.
 */
static double density_upper_bound(density_function_t f,
                                  const prior_lookup_t* g_mode,
                                  const interval_bound_t* bound,
                                  idx_t k_min, idx_t k_max)
{
    if (bound->end_min > bound->start_max) {
        k_min = idxmax(k_min, bound->end_min - bound->start_max + 1);
    }

    k_max = idxmin(k_max, bound->end_max - bound->start_min + 1);

    return f(bound->x_max, k_min) + prior_lookup_eval(g_mode, k_max);
}


/* Parameters passed to each peakolator thread. */
typedef struct peakolator_ctx_t_
{
    /* Data. */
    const vector_t* vec;

    /* How many positions have been searched. */
    pthread_mutex_t done_mutex;
    size_t done;

    /* Intervals left to search. */
    interval_stack_t* in;

    /* High-density intervals found. */
    interval_stack_t* out;

    /* Allowable lengths for highest-density intervals. */
    idx_t k_min, k_max;

    /* Minimum density for a reported high-denisty region. */
    double min_density;

    /* Density fuction. */
    density_function_t f;

    /* Lookup tables for the prior over interval length. */
    prior_lookup_t* g_lookup;
    prior_lookup_t* g_mode_lookup;

} peakolator_ctx_t;


/* Perform a brute-force search for the highest-density interval. */
static void peakolator_brute(peakolator_ctx_t* ctx,
                             const interval_bound_t* bound,
                             interval_t* best)
{
    vector_sum_t a, b;
    idx_t i, j;

    j = bound->start_min + ctx->k_min - 1;
    if (j < bound->end_min) j = bound->end_min;

    vector_sum_set(&a, ctx->vec, bound->start_min, j);

    for (i = bound->start_min; i <= bound->start_max; ++i) {
        j = i + ctx->k_min - 1;
        if (j < bound->end_min) j = bound->end_min;
        vector_sum_update_end(&a, j);
        /*vector_sum_update_start(&a, i);*/
        if (i > bound->start_min) vector_sum_inc_start(&a);
        memcpy(&b, &a, sizeof(vector_sum_t));

        idx_t j_max = i + ctx->k_max - 1;
        if (j_max > bound->end_max) j_max = bound->end_max;

        /* TODO: Why does this actually occur? */
        if (j > j_max) continue;

        while (true) {
            /*vector_sum_update_end(&b, j);*/
            idx_t k = j - i + 1;
            double density = ctx->f(b.sum, k) + prior_lookup_eval(ctx->g_lookup, k);

            if (density > best->density) {
                best->start = i;
                best->end = j;
                best->density = density;
                best->x = b.sum;
            }

            if (j == j_max) break;
            vector_sum_inc_end(&b);
            ++j;
        }
    }
}


/* Use brute-force search when there are fewer than this many intervals. */
static const uint64_t brute_cutoff = 10000;


/* A single peakolator thread. */
static void* peakolator_thread(void* ctx_)
{
    peakolator_ctx_t* ctx = (peakolator_ctx_t*) ctx_;
    pqueue_t* bounds = pqueue_create();

    interval_t interval;
    interval_bound_t bound, subbound;

    while (interval_stack_pop(ctx->in, &interval)) {
        idx_t k = interval.end - interval.start + 1;
        if (k < ctx->k_min) {
            pthread_mutex_lock(&ctx->done_mutex);
            ctx->done += k;
            pthread_mutex_unlock(&ctx->done_mutex);
            interval_stack_finish_one(ctx->in);
            continue;
        }

        bound.start_min = bound.end_min = interval.start;
        bound.start_max = bound.end_max = interval.end;
        bound.x_max = vector_sum(ctx->vec, interval.start, interval.end);
        bound.density_max = density_upper_bound(ctx->f, ctx->g_mode_lookup,
                                                &bound, ctx->k_min, ctx->k_max);
        if (bound.density_max <= ctx->min_density) {
            pthread_mutex_lock(&ctx->done_mutex);
            ctx->done += k;
            pthread_mutex_unlock(&ctx->done_mutex);
            interval_stack_finish_one(ctx->in);
            continue;
        }

        pqueue_enqueue(bounds, &bound);

        interval_t best;
        best.density = ctx->min_density;

        while (pqueue_dequeue(bounds, &bound)) {
            /* Throw out subsets of the search space that could not possibly
             * hold the high density interval. */
            if (bound.density_max <= best.density) {
                break;
            }

            /*if (bound.density_max <= best.density) continue;*/

            uint64_t count = interval_bound_count(&bound,
                                                  ctx->k_min, ctx->k_max);
            if (count == 0)  continue;

            /* Resort to a brute force search when the bound contains relatively
             * few intervals. */
            else if (count < brute_cutoff) {
                peakolator_brute(ctx, &bound, &best);
                continue;
            }

            /* The search space in partitioned into three or four subsets by
             * bisecting the the start and end bounds and taking the cartesian
             * product.
             *
             * Schematically, this looks like so:
             *
             *                Start Bound      End Bound
             * Partition 1.   XXXXXX------     ------XXXXXX
             * Partition 2.   XXXXXX------     XXXXXX------
             * Partition 3.   ------XXXXXX     ------XXXXXX
             * Partition 4.   ------XXXXXX     XXXXXX------
             *
             * The fourth partition is ill-defined is not included when the
             * start bound and end bound are the same (hence there being
             * sometimes three partitions.).
             * */

            /*bool equal_bound = bound.start_min == bound.end_min &&*/
                               /*bound.start_max == bound.end_max;*/

            idx_t start_bound_len = bound.start_max - bound.start_min + 1;
            idx_t end_bound_len   = bound.end_max   - bound.end_min + 1;
            idx_t start_mid = bound.start_min + start_bound_len / 2;
            idx_t end_mid = bound.end_min + end_bound_len / 2;

            val_t x_start0 = vector_sum(ctx->vec, bound.start_min, start_mid);
            val_t x_end1   = vector_sum(ctx->vec, end_mid + 1, bound.end_max);

            /* Partition 1 */
            subbound.start_min   = bound.start_min;
            subbound.start_max   = start_mid;
            subbound.end_min     = end_mid + 1;
            subbound.end_max     = bound.end_max;
            subbound.x_max       = bound.x_max;
            subbound.density_max = density_upper_bound(
                    ctx->f, ctx->g_mode_lookup,
                    &subbound, ctx->k_min, ctx->k_max);
            if (interval_bound_min_len(&subbound) <= ctx->k_max &&
                subbound.density_max > ctx->min_density) {
                pqueue_enqueue(bounds, &subbound);
            }

            /* Partition 2 */
            subbound.start_min   = bound.start_min;
            subbound.start_max   = start_mid;
            subbound.end_min     = bound.end_min;
            subbound.end_max     = end_mid;
            subbound.x_max       = bound.x_max - x_end1;
            subbound.density_max = density_upper_bound(
                    ctx->f, ctx->g_mode_lookup,
                    &subbound, ctx->k_min, ctx->k_max);
            if (interval_bound_min_len(&subbound) <= ctx->k_max &&
                subbound.density_max > ctx->min_density) {
                pqueue_enqueue(bounds, &subbound);
            }

            /* Partition 3 */
            subbound.start_min   = start_mid + 1;
            subbound.start_max   = bound.start_max;
            subbound.end_min     = end_mid + 1;
            subbound.end_max     = bound.end_max;
            subbound.x_max       = bound.x_max - x_start0;
            subbound.density_max = density_upper_bound(
                    ctx->f, ctx->g_mode_lookup,
                    &subbound, ctx->k_min, ctx->k_max);
            if (interval_bound_min_len(&subbound) <= ctx->k_max &&
                subbound.density_max > ctx->min_density) {
                pqueue_enqueue(bounds, &subbound);
            }

            /* Partition 4 */
            if (bound.end_min >= bound.start_max) {
                subbound.start_min   = start_mid + 1;
                subbound.start_max   = bound.start_max;
                subbound.end_min     = bound.end_min;
                subbound.end_max     = end_mid;
                subbound.x_max       = bound.x_max - x_start0 - x_end1;
                subbound.density_max = density_upper_bound(
                        ctx->f, ctx->g_mode_lookup,
                        &subbound, ctx->k_min, ctx->k_max);
                if (interval_bound_min_len(&subbound) <= ctx->k_max &&
                    subbound.density_max > ctx->min_density) {
                    pqueue_enqueue(bounds, &subbound);
                }
            }
        }

        pqueue_clear(bounds);

        pthread_mutex_lock(&ctx->done_mutex);

        /* Find anything good? */
        if (best.density > ctx->min_density) {
            interval_stack_push(ctx->out, &best);
            ctx->done += best.end - best.start + 1;

            idx_t start = interval.start;
            idx_t end   = interval.end;

            if (best.start - start >= ctx->k_min) {
                interval.start = start;
                interval.end   = best.start - 1;
                interval_stack_push(ctx->in, &interval);
            }
            else ctx->done += best.start - start;

            if (end - best.end >= ctx->k_min) {
                interval.start = best.end + 1;
                interval.end   = end;
                interval_stack_push(ctx->in, &interval);
            }
            else ctx->done += end - best.end;
        }
        else ctx->done += interval.end - interval.start + 1;

        pthread_mutex_unlock(&ctx->done_mutex);
        interval_stack_finish_one(ctx->in);
    }

    pqueue_free(bounds);
    return NULL;
}


/* Print a progress indicator in a seperate thread.
 *
 * Args:
 *   ctx_: A pointer to the peakolator_cxt_t that was passed to the peakolator
 *         threads.
 * */
void* peakolator_progress(void* ctx_)
{
    peakolator_ctx_t* ctx = (peakolator_ctx_t*) ctx_;

    size_t barlen = 60;

    fprintf(stderr, "  |\033[s");

    while (!interval_stack_finished(ctx->in)) {
        double p = (double) ctx->done / (double) ctx->vec->n;
        size_t i;
        fprintf(stderr, "\033[u\033[K");

        size_t steps = (size_t) round(p * (double) barlen);
        for (i = 0; i < steps; ++i) fputc('.', stderr);
        for (i = 0; i < barlen - steps; ++i) fputc(' ', stderr);
        fprintf(stderr, "| %0.2f%%", 100.0 * p);

        usleep(100000);
    }

    return NULL;
}


/* The peakolator algorithm.
 *
 * The core the peakolator algorithm is actually embarrassingly simple and
 * included mostly in `peakolator_thread`. All this code is a lot of fuss around
 * making it run efficiently in parallel. This function launches `num_threads`
 * threads which all share one priority queue. The priority queue holds interval
 * bounds (subsets of the search space) in order of maximum density.
 *
 * Each thread dequeues a bound, and either does a brute force search (if the
 * bound is small) or further partitions the bound into subbounds and enqueues
 * these subounds back on the priority queue.
 */
size_t peakolate(const vector_t* vec,
                 density_function_t f,
                 prior_function_t g,
                 idx_t k_min,
                 idx_t k_max,
                 double min_density,
                 unsigned int num_threads,
                 bool show_progress,
                 interval_t** out)
{
    if (num_threads == 0) num_threads = ncpus();

    prior_lookup_t* g_lookup = memoize_prior(g, k_min, k_max);
    prior_lookup_t* g_mode_lookup = memoize_prior_mode(g, k_min, k_max);

    peakolator_ctx_t ctx;
    ctx.vec = vec;
    pthread_mutex_init_or_die(&ctx.done_mutex, NULL);
    ctx.done = 0;
    ctx.in  = interval_stack_create();
    ctx.out = interval_stack_create();
    ctx.k_min = k_min;
    ctx.k_max = k_max;
    ctx.min_density = min_density;
    ctx.f = f;
    ctx.g_lookup = g_lookup;
    ctx.g_mode_lookup = g_mode_lookup;

    /* Initial interval to search. */
    interval_t interval = {0, vector_len(vec) - 1, -INFINITY, 0};
    interval_stack_push(ctx.in, &interval);

    pthread_t* threads = malloc_or_die(num_threads * sizeof(pthread_t));

    size_t i;
    for (i = 0; i < num_threads; ++i) {
        pthread_create(&threads[i], NULL, peakolator_thread, &ctx);
    }

    if (show_progress) {
        pthread_t progress_thread;
        pthread_create(&progress_thread, NULL, peakolator_progress, &ctx);
        pthread_join(progress_thread, NULL);
        fputs("\033[G\033[K", stderr); /* Clear progress bar. */
    }

    for (i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    *out = malloc_or_die(ctx.out->n * sizeof(interval_t));
    memcpy(*out, ctx.out->xs, ctx.out->n * sizeof(interval_t));
    /*sort_intervals_des_density(*out, ctx.out->n);*/
    sort_intervals_asc_pos(*out, ctx.out->n);
    size_t out_len = ctx.out->n;
    pthread_mutex_destroy(&ctx.done_mutex);
    interval_stack_free(ctx.in);
    interval_stack_free(ctx.out);
    prior_lookup_free(g_mode_lookup);
    prior_lookup_free(g_lookup);
    free(threads);

    return out_len;
}


