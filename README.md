
# Peakolator

The world's fastest genomic segmentation algorithm.


# Problem

Peakolator solves the following problem:

Given a vector `(x_1, x_2, ..., x_n)`, find the highest density interval `[i,
j]`, where density is defined by some function `f(x, k)`, where `x = x_i + ... +
x_j` and `k = j - i + 1`. Typically `n` is quite large, say 100 million or so.

The density function must comply with some simple properties and in exchange,
Peakolator can typically arrive the solution to this problem quite quickly.


## Density function

The algorithm centers around a density function (my terminology). A density
function is any function `f(k, x)` that obeys the following two properties.
1. `f(x', k) >= f(x, k)` for all `x' >= x`.
2. `f(x, k') <= f(x, k)` for all `k' >= k`.

The first property says that density is non-decreasing with mass (`x`), and the
second that it is non-increasing with volume (`k`). This should hopefully fit
with one's intuitive notion of "density".

There is an obvious density function (namely, `f(x, k) = x/k`) that obeys these
two properties, but many other functions do as well.


## Scan Statistics

One example is the traditional frequentist scan statistics methodologies.
Suppose we have a vector `(x_1, ..., x_n)` and a null hypothesis that `x_i ~
Poisson(a)`, and thus `x_i + ... + i_j ~ Poisson((j - i + 1) * a)`. We wish to
find intervals in the vector in which the null hypothesis has very
low-likelihood.

This sort of approach is used very frequently with "sliding window" algorithms,
but sliding window algorithms involve choosing one or more fixed window lengths.
Peakolator can search for significant windows of all length.  We can set `f(x,
k) = 1 - P(x | k * a)` where `P` is the Poisson CDF function.  This is a proper
density function and can be used with Peakolator!


## Extensions

A wide variety of useful functions are density functions, but it is also a
limiting notion. Often we have some idea of the length of the regions we are
search for. Peakolator thus supports specifying an arbitrary prior distribution
over interval lengths.


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


