#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Some functions needed by multiple example programs. */


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

/* Simple encapsulated string. */
typedef struct
{
    unsigned char* s;    /* null-terminated string */
    uint32_t       n;    /* length of s */
    uint32_t       size; /* bytes allocated for s */
} str_t;


/* Initialize an empty string. */
void str_init(str_t* str)
{
    str->n    = 0;
    str->size = 0;
    str->s    = NULL;
}


/* Free a string. */
void str_free(str_t* str)
{
    if (str) free(str->s);
}


/* Reserve space for `size` more characters. */
void str_reserve_extra(str_t* str, size_t size)
{
    if (str->n + size > str->size) {
        str->size = str->n + size;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}


/* Append a char to a string. */
void str_append_char(str_t* a, char c)
{
    str_reserve_extra(a, 1);
    a->s[a->n++] = c;
}


/* A sequence with an associated name. */
typedef struct
{
    char* name;
    char* seq;
} namedseq_t;


/* Does a character represent a nucleotide? */
static bool is_nt_char(char c)
{
    return c == 'a' || c == 'A' ||
           c == 'c' || c == 'C' ||
           c == 'g' || c == 'G' ||
           c == 't' || c == 'T' ||
           c == 'n' || c == 'N';
}


/* Read the next sequence in an open FASTA files.
 *
 * Args:
 *   file: An open file with a position before the next entry to be read.
 *   seq: A pointer to pointer which will be set to the sequence.
 *
 * Returns:
 *   true if a new sequence was read.
 *
 */
size_t read_fasta(FILE* file, namedseq_t** out)
{
    size_t size = 16;
    size_t n = 0;
    namedseq_t* seqs = malloc(size * sizeof(namedseq_t));

    const size_t bufsize = 1024;
    char* buf = malloc_or_die(bufsize); buf[0] = '\0';
    char* next = buf;

    str_t name;
    str_t seq;

    str_init(&name);
    str_init(&seq);

    /* The three parser states are:
         0 : reading seqname
         1 : reading sequence
         2 : reading sequence (line beginning)
    */
    int state = 2;

    while (true) {
        /* end of buffer */
        if (*next == '\0') {
            if (fgets(buf, bufsize, file) == NULL) {
                /* end of file */
                break;
            }
            else {
                next = buf;
                continue;
            }
        }
        else if (state == 0) {
            if (*next == '\n') {
                str_append_char(&name, '\0');
                char* c = strchr((char*) name.s, ' ');
                if (c != NULL) {
                    *c = '\0';
                    name.n = c - (char*) name.s;
                }

                fprintf(stderr, "\treading %s...\n", name.s);

                state = 2;

                if (n >= size) {
                    size *= 2;
                    seqs = realloc_or_die(seqs, size * sizeof(namedseq_t));
                }

                seqs[n].name = malloc(name.n + 1);
                memcpy(seqs[n].name, name.s, name.n);
                seqs[n].name[name.n] = '\0';

                ++n;
            }
            else {
                str_append_char(&name, *next);
            }

            ++next;
        }
        else if (state == 1 || state == 2) {
            if (*next == '\n') {
                state = 2;
            }
            else if (is_nt_char(*next)) {
                str_append_char(&seq, *next);
                state = 1;
            }
            else if (state == 2 && *next == '>') {
                if (n > 0) {
                    seqs[n-1].seq = malloc(seq.n + 1);
                    memcpy(seqs[n-1].seq, seq.s, seq.n);
                    seqs[n-1].seq[seq.n] = '\0';
                }

                name.n = 0;
                seq.n = 0;
                state = 0;
            }
            else {
                fprintf(stderr, "Error parsing FASTA file: "
                                "unexpected character '%c'.", *next);
            }

            ++next;
        }
    }

    if (n > 0) {
        seqs[n-1].seq = malloc(seq.n + 1);
        memcpy(seqs[n-1].seq, seq.s, seq.n);
        seqs[n-1].seq[seq.n] = '\0';
    }

    *out = seqs;
    return n;
}


