#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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


/* Compare two interval bounds;
 *
 * Args:
 *   a: An interval bound.
 *   b: Another interval bound.
 *
 * Returns:
 *   Returns 0 if the bounds are equal and non-zero otherwise.
 *   (Hint: this is just a wrapper around memcmp.)
 */
int interval_bound_cmp(const interval_bound_t* a, const interval_bound_t* b);


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
    memcpy(dest, src, sizeof(interval_bound_t));
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
    memcpy(dest, src, sizeof(interval_t));
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


/* A stack of intervals. */
typedef struct interval_stack_t_
{
    /* Items in the stack. */
    interval_t* xs;

    /* Number of items. */
    size_t n;

    /* Size reserved. */
    size_t size;

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
    s->size = 1024;
    s->xs = malloc_or_die(s->size * sizeof(interval_t));
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
    if (s->n == s->size) interval_stack_expand(s);
    interval_copy(&s->xs[s->n++], interval);
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
    if (s->n == 0) return false;
    interval_copy(interval, &s->xs[--s->n]);
    return true;
}


/* A concurrent priority queue for genomic intervals.
 *
 * This is just a coarsely locked binary max-heap, nothing fancy.
 */
typedef struct pqueue_t_
{
    /* Coarse lock. */
    pthread_mutex_t mutex;

    /* Signal waiters that an item is available. */
    pthread_cond_t cond;

    /* Number of items in the heap. */
    size_t n;

    /* Size reserved. */
   size_t size;

    /* Items. */
    interval_bound_t* xs;

    /* Number of pending items. */
    size_t pending;
} pqueue_t;


/* Create a new priority queue. */
static pqueue_t* pqueue_create()
{
    pqueue_t* q = malloc_or_die(sizeof(pqueue_t));
    q->n = 0;
    q->size = 1024;
    q->xs = malloc_or_die(q->size * sizeof(interval_bound_t));
    q->pending = 0;
    pthread_mutex_init_or_die(&q->mutex, NULL);
    pthread_cond_init_or_die(&q->cond, NULL);
    return q;
}


/* Free an priority queue created with pqueue_create. */
static void pqueue_free(pqueue_t* q)
{
    pthread_mutex_destroy(&q->mutex);
    pthread_cond_destroy(&q->cond);
    free(q->xs);
    free(q);
}


/* Double the size reserved for a heap. */
static void pqueue_expand(pqueue_t* q)
{
    q->size *= 2;
    q->xs = realloc_or_die(q->xs, q->size * sizeof(interval_bound_t));
}


/* Binary heap index function. */
static inline size_t pqueue_parent_idx (size_t i) { return (i - 1) / 2; }
static inline size_t pqueue_left_idx   (size_t i) { return 2 * i + 1; }
static inline size_t pqueue_right_idx  (size_t i) { return 2 * i + 2; }


/* Insert an item into the priority queue.
 *
 * Args:
 *   q: A pqueue.
 *   interval: Interval to enqueue.
 */
static void pqueue_enqueue(pqueue_t* q, const interval_bound_t* bound)
{
    pthread_mutex_lock(&q->mutex);

    /* insert */
    if (q->n == q->size) pqueue_expand(q);
    size_t j, i = q->n++;
    interval_bound_copy(&q->xs[i], bound);

    /* percolate up */
    while (i > 0) {
        j = pqueue_parent_idx(i);
        if (q->xs[j].density_max < q->xs[i].density_max) {
            interval_bound_swap(&q->xs[i], &q->xs[j]);
            i = j;
        }
        else break;
    }

    pthread_mutex_unlock(&q->mutex);
    pthread_cond_signal(&q->cond);
}


/* Pop the item with the largest density from the queue.
 *
 * If the queue is empty, this function will block until an item becomes
 * available, unless pqueue_finish has been called, in which case it will return
 * immediately.
 *
 * Args:
 *   q: A pqueue.
 *   interval: Location to copy popped interval.
 *
 * Returns:
 *   false if the queue is finished and empty.
 */
static bool pqueue_dequeue(pqueue_t* q, interval_bound_t* bound)
{
    pthread_mutex_lock(&q->mutex);

    while (q->n == 0 && q->pending > 0) {
        pthread_cond_wait(&q->cond, &q->mutex);
    }

    if (q->n == 0 && q->pending == 0) {
        pthread_mutex_unlock(&q->mutex);
        return false;
    }

    interval_bound_copy(bound, &q->xs[0]);

    /* replace head */
    interval_bound_copy(&q->xs[0], &q->xs[--q->n]);

    /* percolate down */
    size_t l, r, j, i = 0;
    while (true) {
        l = pqueue_left_idx(i);
        r = pqueue_right_idx(i);

        if (l < q->n) {
            j = r < q->n && q->xs[r].density_max > q->xs[l].density_max ? r : l;

            if (q->xs[j].density_max > q->xs[i].density_max) {
                interval_bound_swap(&q->xs[j], &q->xs[i]);
                i = j;
            }
            else break;
        }
        else break;
    }

    ++q->pending;
    pthread_mutex_unlock(&q->mutex);
    return true;
}


/* Decrement the number of pending items.
 *
 * If empty the queue blocks on dequeue unless the number of pending items is
 * zero. That number is incremented on dequeue. The caller of dequeue has to call
 * this function at a point where it will not enqueue any more without first
 * dequeuing.
 *
 * Args:
 *  q: A pqueue.
 */
static void pqueue_finish_one(pqueue_t* q)
{
    pthread_mutex_lock(&q->mutex);
    if (q->pending > 0) --q->pending;
    pthread_mutex_unlock(&q->mutex);

    // In case anyone is waiting on an item that will never arrive.
    pthread_cond_broadcast(&q->cond);
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
    prior_lookup_t* lookup =
        malloc_or_die(sizeof(prior_lookup_t));
    lookup->xs = malloc_or_die((max - min + 1) * sizeof(double));
    lookup->min = min;
    lookup->max = max;

    idx_t i;
    for (i = min; i <= max; ++i) {
        lookup->xs[i] = g(i);
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
    prior_lookup_t* lookup =
        malloc_or_die(sizeof(prior_lookup_t));
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
                                  idx_t k_min,
                                  const prior_lookup_t* g_mode,
                                  val_t x, idx_t k)
{
    return f(x, k_min) + prior_lookup_eval(g_mode, k);
}


/* Parameters passed to each peakolator thread. */
typedef struct peakolator_ctx_t_
{
    /* Data. */
    const vector_t* vec;

    /* Fragments left to search. */
    pqueue_t* in;

    /* Current highest density interval. Density is -Inf if none has been found
     * yet. */
    interval_t best;

    /* Mutex for atomic access to "best". */
    pthread_mutex_t best_mutex;

    /* Allowable lengths for highest-density intervals. */
    idx_t k_min, k_max;

    /* Density fuction. */
    density_function_t f;

    /* Lookup tables for the prior over interval length. */
    prior_lookup_t* g_lookup;
    prior_lookup_t* g_mode_lookup;

} peakolator_ctx_t;


/* Perform a brute-force search for the highest-density interval. */
static void peakolator_brute(peakolator_ctx_t* ctx,
                             const interval_bound_t* bound)
{
    double best_density = ctx->best.density;

    idx_t i, j;
    for (i = bound->start_min; i <= bound->start_max; ++i) {
        j = i + ctx->k_min - 1;
        if (j < bound->end_min) j = bound->end_min;

        idx_t j_max = i + ctx->k_max - 1;
        if (j_max > bound->end_max) j_max = bound->end_max;

        for (; j <= j_max; ++j) {
            /* TODO:
             * We could be more clever with counting, since we already
             * know the sum of [i, j - 1].
             *
             * Actually, we really need to be. About 90% of the runtime is the
             * line below. */
            val_t x = vector_sum(ctx->vec, i, j);
            idx_t k = j - i + 1;
            double density = ctx->f(x, k) + prior_lookup_eval(ctx->g_lookup, k);

            if (density > best_density) {
                pthread_mutex_lock(&ctx->best_mutex);
                if (density > ctx->best.density) {
                    ctx->best.start = i;
                    ctx->best.end = j;
                    ctx->best.density = density;
                }
                best_density = ctx->best.density;
                pthread_mutex_unlock(&ctx->best_mutex);
            }
        }
    }
}


/* Use brute-force search when there are fewer than this many intervals. */
static const uint64_t brute_cutoff = 1000;


/* A single peakolator thread. */
static void* peakolator_thread(void* ctx_)
{
    peakolator_ctx_t* ctx= (peakolator_ctx_t*) ctx_;
    interval_bound_t bound, subbound;

    while (pqueue_dequeue(ctx->in, &bound)) {
        uint64_t count = interval_bound_count(&bound, ctx->k_min, ctx->k_max);

        if (count == 0) {
            pqueue_finish_one(ctx->in);
            continue;
        }
        else if (count < brute_cutoff) {
            peakolator_brute(ctx, &bound);
            pqueue_finish_one(ctx->in);
            continue;
        }

        /* The search space in partitioned into three or four subsets by
         * bisecting the the start and end bounds and taking the cartesian
         * product.
         *
         * Schematically, this looks like so:
         *
         *                Start Bound      End Bound
         * Partition 1.   ######------     ------######
         * Partition 2.   ######------     ######------
         * Partition 3.   ------######     ------######
         * Partition 4.   ------######     ######------
         *
         * The fourth partition is ill-defined is not included when the start
         * bound and end bound are the same (hence there being sometimes three
         * partitions.).
         * */

        bool equal_bound = bound.start_min == bound.end_min &&
                           bound.start_max == bound.end_max;

        idx_t start_bound_len = bound.start_max - bound.start_min + 1;
        idx_t end_bound_len   = bound.end_max   - bound.end_min + 1;
        idx_t start_mid = bound.start_min + start_bound_len / 2;
        idx_t end_mid = bound.end_min + end_bound_len / 2;

        val_t x_start0 = vector_sum(ctx->vec, bound.start_min, start_mid);
        val_t x_start1 = vector_sum(ctx->vec, start_mid + 1, bound.start_max);
        val_t x_end0   = vector_sum(ctx->vec, bound.end_min, end_mid);
        val_t x_end1   = vector_sum(ctx->vec, end_mid + 1, bound.end_max);

        /* Partition 1 */
        subbound.start_min   = bound.start_min;
        subbound.start_max   = start_mid;
        subbound.end_min     = end_mid + 1;
        subbound.end_max     = bound.end_max;
        subbound.x_max       = x_start0 + x_end1;
        subbound.density_max = density_upper_bound(
                ctx->f, ctx->k_min, ctx->g_mode_lookup,
                subbound.x_max, subbound.end_max - subbound.start_min + 1);
        pqueue_enqueue(ctx->in, &subbound);

        /* Partition 2 */
        subbound.start_min   = bound.start_min;
        subbound.start_max   = start_mid;
        subbound.end_min     = bound.end_min;
        subbound.end_max     = end_mid;
        subbound.x_max       = equal_bound ? x_start0 : x_start0 + x_end0;
        subbound.density_max = density_upper_bound(
                ctx->f, ctx->k_min, ctx->g_mode_lookup,
                subbound.x_max, subbound.end_max - subbound.start_min + 1);
        pqueue_enqueue(ctx->in, &subbound);

        /* Partition 3 */
        subbound.start_min   = start_mid + 1;
        subbound.start_max   = bound.start_max;
        subbound.end_min     = end_mid + 1;
        subbound.end_max     = bound.end_max;
        subbound.x_max       = equal_bound ? x_start1 : x_start1 + x_end1;
        subbound.density_max = density_upper_bound(
                ctx->f, ctx->k_min, ctx->g_mode_lookup,
                subbound.x_max, subbound.end_max - subbound.start_min + 1);
        pqueue_enqueue(ctx->in, &subbound);

        /* Partition 4 */
        if (bound.end_min >= bound.start_max) {
            subbound.start_min   = start_mid + 1;
            subbound.start_max   = bound.start_max;
            subbound.end_min     = bound.end_min;
            subbound.end_max     = end_mid;
            subbound.x_max       = x_start1 + x_end0;
            subbound.density_max = density_upper_bound(
                    ctx->f, ctx->k_min, ctx->g_mode_lookup,
                    subbound.x_max, subbound.end_max - subbound.start_min + 1);
            pqueue_enqueue(ctx->in, &subbound);
        }

        pqueue_finish_one(ctx->in);
    }

    return NULL;
}

/* TODO: Describe what it is I'm doing here.
 */
void peakolate(const vector_t* vec,
               density_function_t f,
               prior_function_t g,
               idx_t k_min,
               idx_t k_max,
               int num_threads)
{
    /* TODO: if num_threads <= 0 set it to be equal to the number of cores. */

    pthread_t* threads = malloc_or_die(num_threads * sizeof(pthread_t));

    prior_lookup_t* g_lookup = memoize_prior(g, k_min, k_max);
    prior_lookup_t* g_mode_lookup = memoize_prior_mode(g, k_min, k_max);

    peakolator_ctx_t ctx;
    ctx.vec = vec;
    ctx.in = pqueue_create();
    pthread_mutex_init_or_die(&ctx.best_mutex, NULL);
    ctx.k_min = k_min;
    ctx.k_max = k_max;
    ctx.f = f;
    ctx.g_lookup = g_lookup;
    ctx.g_mode_lookup = g_mode_lookup;

    /* Queue of intervals left to search. */
    interval_stack_t* s = interval_stack_create();

    interval_t interval;
    interval.start = 0;
    interval.end = vector_len(vec) - 1;
    interval_stack_push(s, &interval);

    interval_bound_t bound;
    idx_t k;

    while (interval_stack_pop(s, &interval)) {
        k = interval.end - interval.start + 1;
        if (k < k_max) continue;

        /* Enqueue the initial interval bound. */
        bound.start_min = bound.end_min = interval.start;
        bound.start_max = bound.end_max = interval.end;
        bound.x_max = vector_sum(vec, interval.start, interval.end);
        bound.density_max = density_upper_bound(f, k_min, g_mode_lookup,
                                                bound.x_max, k);
        pqueue_enqueue(ctx.in, &bound);

        ctx.best.density = -INFINITY;

        /* Start up a bunch of search threads. */
        int i;
        for (i = 0; i < num_threads; ++i) {
            pthread_create(&threads[i], NULL, peakolator_thread, &ctx);
        }

        for (i = 0; i < num_threads; ++i) {
            pthread_join(threads[i], NULL);
        }

        /* Find anything good? */
        if (ctx.best.density != -INFINITY) {
            /* TODO: Do something with ctx.best when we figure
             * out what this function returs.. */

            if (ctx.best.start > interval.start) {
                interval.end = ctx.best.start - 1;
                interval_stack_push(s, &interval);
            }

            if (ctx.best.end < interval.end) {
                interval.start = ctx.best.end + 1;
                interval_stack_push(s, &interval);
            }
        }
    }

    interval_stack_free(s);
    pqueue_free(ctx.in);
    pthread_mutex_destroy(&ctx.best_mutex);
    prior_lookup_free(g_mode_lookup);
    prior_lookup_free(g_lookup);
    free(threads);
}


