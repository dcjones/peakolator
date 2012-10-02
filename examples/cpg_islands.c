
#include <ctype.h>
#include <getopt.h>
#include <math.h>

#include "common.h"
#include "../src/peakolator.h"


static idx_t min_length  = 500;
static idx_t max_length  = 10000;
static double mean_length = 1000;
static double std_length  = 500;
static double fore_log_p, fore_log_q, back_log_p, back_log_q;


static double logaddexp(double x, double y)
{
    double u = x - y;
    if (u > 0.0) {
        return x + log1p(exp(-u));
    }else if (u <= 0.0) {
        return y + log1p(exp(u));
    }else  {
        return x + y;
    }
}


double f(val_t x_, idx_t k_)
{
/* This objective function is fucked. It's very easy to arrive at a posterion
 * probability of 1.0 */

    double x = (double) x_;
    double k = (double) k_;

    if (x > k) x = k;

    return x / k;
#if 0

    /* These are binomial distribution propabilities, but witohut the binomial
     * coefficient term, since that drops out of the posterior probality. */
    double u = x * fore_log_p + (k - x) * fore_log_q;
    double v = x * back_log_p + (k - x) * back_log_q;

    /* TODO: prior */
    if (u - logaddexp(u, v) == 0) {
        fprintf(stderr, "HERE\n");
    }

    return u - logaddexp(u, v);
#endif
}


double g(__attribute__((unused)) idx_t k)
{
    /* TODO */
    return 0;
}


void print_usage(FILE* file)
{
    fprintf(file,
"Usage: cpg_islands [OPTIONS]\n"
"Find CpG islands in a set of nucleotide sequences.\n"
"Sequences are read from standard input in FASTA format, and predicted CpG\n"
"islands are outputted in BED format.\n"
"\n"
"Options:\n");
    /* TODO */
}


int main(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"min-length",  required_argument, NULL, 'm'},
        {"max-length",  required_argument, NULL, 'M'},
        {"mean-length", required_argument, NULL, 'u'},
        {"std-length",  required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}
    };


    int opt, opt_idx;
    while (true) {
        opt = getopt_long(argc, argv, "h", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'm':
                min_length = atoi(optarg);
                break;

            case 'M':
                max_length = atoi(optarg);
                break;

            case 'u':
                mean_length = atof(optarg);
                break;

            case 's':
                std_length = atof(optarg);
                break;

            case 'h':
                print_usage(stdout);
                return EXIT_SUCCESS;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }

    }

    namedseq_t* seqs;
    size_t num_seqs = read_fasta(stdin, &seqs);
    size_t i, j;

    /* Determine background CG rate. */
    uint64_t cg_count = 0;
    uint64_t count = 0;
    for (i = 0; i < num_seqs; ++i) {
        size_t n = strlen(seqs[i].seq);
        for (j = 0; j < n - 1; ++j) {
            if (toupper(seqs[i].seq[j])     == 'C' &&
                toupper(seqs[i].seq[j + 1]) == 'G') {
                cg_count++;
            }
        }
        count += n - 1;
    }
    double cg_pr = (double) cg_count / (double) count;
    fprintf(stderr, "background CG probability: %0.4f\n", cg_pr);

    back_log_p = log(cg_pr);
    back_log_q = log(1.0 - cg_pr);

    cg_pr = fmin(10.0 * cg_pr, 1.0);
    fore_log_p = log(cg_pr);
    fore_log_q = log(1.0 - cg_pr);

    /* Allocate a vector. */
    size_t xs_size = 0;
    for (i = 0; i < num_seqs; ++i) {
        size_t n = strlen(seqs[i].seq);
        if (n > xs_size) xs_size = n;
    }
    double* xs = malloc_or_die(xs_size * sizeof(double));

    /* Process sequences */
    for (i = 0; i < num_seqs; ++i) {
        fprintf(stderr, "searching %s...\n", seqs[i].name);

        size_t n = strlen(seqs[i].seq);
        if (n < 2) continue;

        memset(xs, 0, n * sizeof(double));

        size_t j;
        for (j = 0; j < n - 1; ++j) {
            if (toupper(seqs[i].seq[j])     == 'C' &&
                toupper(seqs[i].seq[j + 1]) == 'G') {
                xs[j] = 1.0;
            }
        }

        vector_t* vec = vector_create(xs, n);

        interval_t* out;
        size_t out_count = peakolate(vec, f, g,
                                     min_length, max_length,
                                     0.10,
                                     /*log(0.99),*/
                                     0, &out);

        for (j = 0; j < out_count; ++j) {
            printf("%s\t%lu\t%lu\t%e\n", seqs[i].name,
                   (unsigned long) out[j].start, (unsigned long) out[j].end + 1,
                   out[j].x);
                   /*out[j].density / M_LN10);*/
        }

        free(out);
        vector_free(vec);
    }

    return EXIT_SUCCESS;
}

