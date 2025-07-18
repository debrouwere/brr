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

### Usage

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

#### Usage with long format imputed data (beta)

The function `brrl` allows for zero-overhead analysis of long format PISA data
in which not only the plausible values are imputed but also the covariates.

```r
fits <- brrl(
  statistic = \(data, weights) lm(math ~ as.factor(cycle), data = data, weights = weights),
  data = assessments_long,
  weights = wts,
  conditions = list(
    i = 1:10,
    weights = c('w_math_student_final', 'w_math_student_r{1:80}')
  )
)

coef(fits)
confint(fits, na_rm = TRUE)
```

Because `weights` are shared by all imputations, they are typically provided
by a data frame with 1/5th or 1/10th the amount of rows as the `data` set.

See `pisa.rx.parquet` for a multiply imputed dataset that is compatible with
`brrl`, or alternatively you can prepare such a dataset yourself using
`dplyr::pivot_longer` and `mice::mice`.

```r
# this is just a sketch, actual conversion into long form may be more involved
# depending on the exact dataset under consideration
outcomes <- assessments |> 
    select(starts_with('pv')) |> 
    pivot_longer(
      cols = everything(), 
      names_pattern = 'pv(\\d+)(math)', 
      names_to = c('i', '.value')
    )

imputations <- assessments |> 
  select(c(escs, fathers_isced, mothers_isced)) |> 
  mice(m = 5) |> 
  complete(action = "long")

data <- bind_cols(outcomes, imputations)

weights <- assessments |> 
  select(starts_with('w_'))
```

In the future, `brrl` will likely become the default analysis interface and the `brr` function
will be deprecated. When `brrl` detects wide format PISA data, it will pivot it
to long format, and therefore it will remain possible to work with the original PISA data 
as provided by the OECD.

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

We often want to poststratify not across the entire dataset but in blocks, e.g.
within a country and an assessment, as shown above. To do so: (1) group by that
block or blocks, (2) poststratify within each block level of each block using
`group_modify` and finally (3) ungroup, e.g. to poststratify separately within
each year of a longitudinal dataset:

```r
belgium03 <- belgium |> filter(assessment == '2003')
belgium |> 
  group_by(year) |>
  group_modify(~ poststratify(.x, poststrata(.x, belgium03, PS_VARS, weights = "w_math_fstuwt")) |>
  ungroup()
```

By default, `poststratify` will accept missingness as a valid stratification
level, so you might get cells like `('Female', '<NA>')`. This is not always a good
idea, both because it it can produce small cells (if there is not much
missingness) and because missingness on a particular factor may have many causes
and does not always represent a homogeneous category. It is often a good idea to
either remove these observations with missingness, or to impute their
stratification factors before poststratification, e.g. with the `mice` package.

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


### Performance

Long format analyses with `brrl` (which will become the future default) are as
fast or faster than alternative software such as EdSurvey, though the older
interfaces such as `brrw` are slower.

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
