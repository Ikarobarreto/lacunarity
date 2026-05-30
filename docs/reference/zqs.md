# Moment generating function of the box masses

Internal helper that computes the \\q\\-th moment \\Z(q) = \sum_m m^q\\
Q(m)\\ of a box-mass frequency distribution, where \\Q(m)\\ is the
relative frequency of boxes carrying mass \\m\\.

## Usage

``` r
zqs(mat, q)
```

## Arguments

- mat:

  a two-column table with the distinct masses `x` and their frequencies
  `freq`.

- q:

  the moment order.

## Value

The \\q\\-th moment of the box-mass distribution.
