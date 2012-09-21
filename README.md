
# Peakolator

The world's fastest genomic segmentation algorithm.


# Problem

Peakolator solves the following problem:

Given a vector `(x_1, x_2, ..., x_n)`, find the highest density interval `[i,
j]`, where density is defined by some function `f(x, k)`, where `x = x_i + ... +
x_j` and `k = j - i + 1`. Typicall `n` is quite large, say 100 million or so.

The density function must comply with some simple properties and in exchange,
Peakolator can typically arrive the solution to this problem.


# Applications

This is a useful thing to be able to do in genomics. First of all, if we can
find the highest density region efficiently, we can perform clustering or
segmentation like so.

1. Find the highest density region.
2. Subtract it from the search space.
3. Repeat.

With the right vector and density function this can be used to solve many
bread-and-butter problems in genomics. E.g. discover CpG islands or highly
conserved regions. It can also be applied to sequencing data to discover novel
transcribed regions from RNA-Seq data, or peaks in ChIP-Seq data.


# Algorithm

Peakolator is a heuristic search algorithm: it has a quadratic worst case
performance, but is almost always extremely fast. It combines two mainstays of
heuristic search: branch-and-bound and best-first.

## Branch and Bound

TODO: Write this.

## Best First

TODO: Write this.


