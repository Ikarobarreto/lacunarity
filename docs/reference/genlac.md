# Generalized lacunarity of a binary series

Computes the generalized (multifractal-like) lacunarity \\\Lambda_q(s)\\
of a binary time series and the spectrum of scaling exponents
\\\gamma(q)\\.

## Usage

``` r
genlac(x)
```

## Arguments

- x:

  a binary vector of 0's and 1's.

## Value

A list with components:

- `s`:

  the dyadic box scales \\s = 2^i\\.

- `q`:

  the moment orders.

- `yq`:

  the generalized scaling exponents \\\gamma(q)\\.

- `Dqs`:

  the matrix of generalized lacunarities \\\Lambda_q(s)\\ (rows index
  `q`, columns index `s`).

## Details

The ordinary lacunarity is extended to an arbitrary moment order \\q\\
by \$\$\Lambda_q(s) = \left\[ \frac{Z(2q,s)}{Z(q,s)^2}
\right\]^{1/q},\$\$ where \\Z(q,s)\\ is the \\q\\-th moment of the
gliding-box mass distribution at scale \\s\\. Large positive \\q\\
emphasises dense boxes and negative \\q\\ emphasises sparse boxes, so
the curve \\q \mapsto \gamma(q)\\, with \\\gamma(q)\\ the slope of
\\\log\_{10}\Lambda_q(s)\\ on \\\log\_{10} s\\, describes how gaps of
different magnitudes scale. Orders \\q\\ range over \\\\-10, \dots, 10\\
\setminus \\0\\\\.

## References

Vernon-Carter, J., Lobato-Calleros, C., Escarela-Perez, R., Rodriguez,
E. and Alvarez-Ramirez, J. (2009). A suggested generalization for the
lacunarity index. *Physica A*, 388(20), 4305-4314.

Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of
random and deterministic fractal sets. *Physical Review A*, 44(6),
3552-3558.

## See also

[`lac`](https://ikarobarreto.github.io/lacunarity/reference/lac.md) for
the ordinary lacunarity index.

## Examples

``` r
x <- rbinom(1000, 1, 0.85)
genlac(x)
#> $s
#> [1]  2  4  8 16 32
#> 
#> $q
#>  [1] -10  -9  -8  -7  -6  -5  -4  -3  -2  -1   1   2   3   4   5   6   7   8   9
#> [20]  10
#> 
#> $yq
#>  [1] -0.129029798 -0.127437583 -0.123564892 -0.117097554 -0.107865779
#>  [6] -0.095671049 -0.079801409 -0.058866781 -0.034187409 -0.013931869
#> [11]  0.008652209  0.013091887  0.015879023  0.017723426  0.019007241
#> [16]  0.019927625  0.020595772  0.021080705  0.021428489  0.021671409
#> 
#> $Dqs
#>            [,1]      [,2]      [,3]      [,4]      [,5]
#>  [1,] 0.8840649 0.6261828 0.6187013 0.6957642 0.8239563
#>  [2,] 0.8724900 0.5962788 0.6181908 0.7092602 0.8401605
#>  [3,] 0.8587241 0.5630036 0.6267368 0.7308015 0.8601418
#>  [4,] 0.8424257 0.5276539 0.6479367 0.7615148 0.8829643
#>  [5,] 0.8235679 0.4944308 0.6853050 0.8010720 0.9069246
#>  [6,] 0.8031368 0.4736261 0.7400805 0.8463629 0.9299439
#>  [7,] 0.7847471 0.4861589 0.8070677 0.8914645 0.9503088
#>  [8,] 0.7778004 0.5674688 0.8737167 0.9305101 0.9672149
#>  [9,] 0.8010085 0.7395506 0.9288656 0.9609317 0.9807400
#> [10,] 0.8709055 0.9085409 0.9700189 0.9834383 0.9914463
#> [11,] 1.0949006 1.0456410 1.0221343 1.0124536 1.0069722
#> [12,] 1.1043947 1.0696495 1.0387693 1.0221187 1.0127857
#> [13,] 1.0966700 1.0817520 1.0513894 1.0298605 1.0177323
#> [14,] 1.0839926 1.0867143 1.0609484 1.0362340 1.0220070
#> [15,] 1.0719088 1.0873179 1.0680987 1.0415987 1.0257414
#> [16,] 1.0618259 1.0853203 1.0733186 1.0461930 1.0290263
#> [17,] 1.0537394 1.0818602 1.0769785 1.0501792 1.0319279
#> [18,] 1.0472946 1.0776715 1.0793756 1.0536711 1.0344969
#> [19,] 1.0421209 1.0732162 1.0807535 1.0567504 1.0367746
#> [20,] 1.0379154 1.0687740 1.0813135 1.0594773 1.0387956
#> 
```
