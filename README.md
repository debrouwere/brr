# brr

Tools for the analysis of data from the Programme for International Student
Assessment (PISA) using balanced repeated replication.

This software is released for educational purposes, not for active use.
It is possible that questions and bug reports may go unanswered.

### Motivation

Although a number of packages exist for the analysis of large-scale educational
assessments, such as `intsvy`, they tend to be limited in the kinds of analyses
they support. This is fine for the needs of educational researcher but
insufficient for statistical research into new or better analysis methods. `brr`
solves this problem by accepting an arbitrary fitter or statistical function to
which it passes a modified formula (for a single plausible value at a time) and
a different set of (replication or final) weights each time, very similar to 
how the bootstrap package `boot` works.

Furthermore, `brr` is written in a succinct and very readable coding style so
that it can be used to explain the technical intricacies of PISA analysis for
those who wish to understand how analyses with plausible values and balanced
repeated replications work under the hood.

## Usage

```r
library('arrow')
library('dplyr')
library('brr')

pisa <- open_dataset('2022.parquet') |>
  select(starts_with('pv') & ends_with('read'), 'country_iso', 'assessment', starts_with('w_read')) |>
  collect()

fits <- brr(pv1read + pv2read + pv3read + pv4read + pv5read + pv6read + pv7read + pv8read + pv9read + pv10read ~ country_iso,
    statistic=lm,
    data=pisa,
    final_weights=matches('w_read_fstuwt'),
    replicate_weights=starts_with('w_read_fstr')
    )

coef(fits)
confint(fits)
confint(fits, extra=TRUE)
```

We use `broom::tidy` to extract standardized output from various kinds of
statistical model.

### Poststratification

This package also implements a poststratification routine similar to the one
used PISA reports (e.g. Chapter 6, 2022 Results, Volume 1) to make trends
comparable over time, by adjusting the demographics of each assessment to match
a target demographical composition.

```r
PS_VARS <- c('age', 'gender', 'immig')

baseline <- pisa |> filter(country_iso == 'FIN', assessment == 2022)

# produce a PISA dataset with poststratified weights
poststratified <- pisa |> 
  group_by(country_iso, assessment) |> 
  group_modify(~ poststratify(.x, poststrata(.x, baseline, PS_VARS, weights = "w_math_fstuwt"), name = "waf")) |> 
  ungroup() |> 
  mutate(across(starts_with("w_"), ~ .x * waf, .names = "poststratified_{col}"))

fits <- brr(pv1math + pv2math + pv3math ~ 1, 
    statistic=weighted_mean_by,
    data=poststratified,
    final_weights=matches('poststratified_w_math_fstuwt'),
    replicate_weights=starts_with('poststratified_w_math_fstr')
    )

confint(fits)
```

### Weighted statistics

Because PISA requires the use of weighted statistics throughout, `brr` comes
bundled with a bunch of weighted statistics with a consistent interface.

* weighted_aggregate
* weighted_mean
* weighted_median
* weighted_quantile
* weighted_sum
* weighted_scale
* weighted_var
* weighted_sd
* weighted_cov
* weighted_cor
* ...

Of particular note is `weighted_aggregate`, which allows for the computation of
univariate statistics for multiple outcomes and multiple subgroups:

```r
df <- tibble(x=1:10, y=11:20, w=1 + 1:10/10, u=rep(1, 10), g=rep(c(1,2), each=5))
weighted_aggregate(x + y ~ g, statistic=weighted_mean, data=df, weights=df$w)
```

For convenience's sake, `weighted_aggregate(..., statistic=weighted_mean)` is
available as `weighted_mean_by`.

A very simple but very fast version that does not support subgroups is available
as `pull_weighted_mean`, very useful if you have already divided the PISA
dataset into subsets (by assessment, by country, ...) elsewhere:

```
results_by_country <- map(countries, function(country) {
  data <- pisa |> filter(country_iso == country)
  brr(pv1scie + pv2scie + pv3scie ~ 1, pull_weighted_mean, ...)
})
```

### Functional helpers

`brr` supports arbitrary fitters and statistical functions, but it does require
them to adhere to a particular interface: `statistic(formula, data, weights,
...)`. If your fitter has a different signature, you can always wrap it in
another function that does, e.g.

```r
wrapped_fitter <- function(formula, data, weights, ...) {
  gnarly(w=weights, d=data, fml=formula)
}
```

The `brr` package also provides a couple of helper functions to cover common
use-cases.

The wrapper `redirect_weights` can be used to rename the `weights` argument
and/or to pass weights as a column name instead of the actual data:

```r
# `brr::brr` passes a labeled weight vector to `statistic`, which we'll emulate here
weights <- label_vector(pisa$w_math_fstuwt, 'w_math_fstuwt')

# adapt a fitter that expects weights to be a column name
lm_by_ref <- function(formula, data, weights) {
  weights <- data[[weights]]
  lm(formula, data, weights=weights)
}

redirect_weights(lm_by_ref, by_reference=TRUE)(pv1read ~ gender, data=pisa, weights=weights)

# adapt a fitter that has a different argument name for weights
lm_wt <- function(formula, data, wt) {
 lm(formula, data, weights=wt)
}

redirect_weights(lm_wt, name='wt')(pv1read ~ gender, data=pisa, weights=weights)
```

The wrapper `formulaic` can add a formula interface to simple statistical
functions:

```r
# this function cannot used in `brr::brr` because it takes data and weight vectors
# instead of a formula and data
weighted_mean(pisa$pv1math, pisa$w_math_fstuwt)

# ... but this one can be used in `brr::brr`!
pull_weighted_mean <- formulaic(weighted_mean)
weighted_mean(pv1math ~ 1, pisa, pisa$w_math_fstuwt)
```

The package includes `pull_weighted_mean` and `pull_weighted_median` out of the
box, but you can use `formulaic` on any weighted statistical function you need.

### Performance

The `brr` package is not particularly fast, but also not particularly slow.

In PISA 2000 up to 2012, proper confidence intervals require 405 balanced
repeated replications: for 5 plausible values and 80 replicate weights plus 1
final weight. Later editions have 10 plausible values, so that's 810 analyses,
which can take a while. However, if only point estimates are required and not
confidence intervals, it is possible to only perform a single analysis for each
plausible value, using the final weights, which is very fast.

For unbiased but approximate point estimates, you can also pick an arbitrary
plausible value and run the analysis outside of `brr`, e.g. using `lm` or
`weighted_mean_by`, taking care to use the final weights.

When setting up your model, it is recommended that you use `brr(r=0, ...)` to
skip the replications, which are only needed for confidence intervals but not
for point estimates. Only when you have a model you are happy with and that runs
without errors, add in the replications.
