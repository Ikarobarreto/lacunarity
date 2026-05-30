# Lacunarity index of a binary series

Computes the gliding-box lacunarity index \\\Lambda(s)\\ of a binary
time series across dyadic scales, together with its scaling exponent.

## Usage

``` r
lac(x)
```

## Arguments

- x:

  a binary vector of 0's and 1's.

## Value

A list with components:

- `y`:

  the lacunarity scaling exponent \\\gamma\\.

- `Ds`:

  the lacunarity \\\Lambda(s)\\ at each scale.

- `s`:

  the dyadic box scales \\s = 2^i\\.

## Details

A box of size \\s\\ is slid one observation at a time along the series
and its mass \\m\\ (the number of ones it covers) is recorded. Writing
\\Z(q,s)\\ for the \\q\\-th moment of the resulting box-mass
distribution, the lacunarity index is \$\$\Lambda(s) =
\frac{Z(2,s)}{Z(1,s)^2} = 1 +
\frac{\mathrm{Var}(m)}{\mathrm{mean}(m)^2},\$\$ so that \\\Lambda(s) \ge
1\\, with equality only for a translationally homogeneous pattern.
Larger values indicate gappier, more heterogeneous textures. The scaling
exponent `y` is the slope of \\\log_2 \Lambda(s)\\ regressed on \\\log_2
s\\. Scales are dyadic, \\s = 2^i\\, and capped by the longest run of
ones.

## References

Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of
random and deterministic fractal sets. *Physical Review A*, 44(6),
3552-3558.

Plotnick, R. E., Gardner, R. H., Hargrove, W. W., Prestegaard, K. and
Perlmutter, M. (1996). Lacunarity analysis: a general technique for the
analysis of spatial patterns. *Physical Review E*, 53(5), 5461-5468.

## See also

[`genlac`](https://ikarobarreto.github.io/lacunarity/reference/genlac.md)
for the generalized lacunarity spectrum.

## Examples

``` r
x <- rbinom(1000, 1, 0.85)
lac(x)
#> $y
#> [1] 0.007806507
#> 
#> $Ds
#> [1] 1.083402 1.042143 1.022033 1.010904 1.005250
#> 
#> $s
#>      [,1]
#> [1,]    2
#> [2,]    4
#> [3,]    8
#> [4,]   16
#> [5,]   32
#> 
```
