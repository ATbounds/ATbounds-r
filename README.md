
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ATbounds

<!-- badges: start -->
<!-- badges: end -->

***ATbounds*** is an R package that provides estimation and inference
methods for bounding average treatment effects (on the treated) that are
valid under an unconfoundedness assumption. The bounds are designed to
be robust in challenging situations, for example, when the the
conditioning variables take on a large number of different values in the
observed sample, or when the overlap condition is violated. This
robustness is achieved by only using limited “pooling” of information
across observations.

For more details, see Lee and Weidner (2021), “Bounding Treatment
Effects by Pooling Limited Information across Observations,” available
at [arXiv:2111.05243 \[econ.EM\]](https://arxiv.org/abs/2111.05243).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") if devtools is not installed
devtools::install_github("sokbae/ATbounds")
```

## Examples

See the vignette for case studies using real datasets.
