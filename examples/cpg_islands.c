
#include <ctype.h>
#include <getopt.h>

#include "common.h"
#include "../src/peakolator.h"


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


    idx_t min_length  = 100;
    idx_t max_length  = 100000;
    double mean_length = 1000;
    double std_length  = 500;

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

        /* TODO: peakolate */

        vector_free(vec);
    }

    return EXIT_SUCCESS;
}

