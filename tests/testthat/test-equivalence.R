library("arrow")
library("dplyr")

PISA_PATH <- Sys.getenv("BRR_PISA_RX_PATH")

EU15_ISO <- c("BEL", "DNK", "DEU", "FIN", "FRA", "GRC", "IRL", "ITA", "LUX", "NLD", "AUT", "PRT", "ESP", "GBR", "SWE")
OECD30_ISO <- c("AUT", "BEL", "CAN", "DNK", "FRA", "DEU", "GRC", "ISL", "IRL", "ITA", "LUX", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "TUR", "GBR", "USA", "JPN", "FIN", "AUS", "NZL", "POL", "HUN", "CZE", "SVK", "KOR", "MEX")

# The equivalence test checks whether numbers computed with `brr` on the student
# data match the per-country summaries published by PISA as XLSX tables attached
# to the publication of the official results.
#
# This test can take a few minutes to run and requires the `pisa.rx.parquet`
# harmonized 2000-2022 dataset, so it is not run as part of the standard test
# suite, but only occasionally to catch regressions after big changes to the
# core `brr::brr` code.
#
# To include this test in a test run from within R, use a `.Renviron` file or:
#
#     Sys.setenv(BRR_TEST_EQUIVALENCE="ON"); devtools::test()
#
# We test to a tolerance of 0.001 (1e-3) at the aggregate level: `expect_equal`
# tests `mean(abs(x - y) < tolerance`; for mean(reading) in 2022 these tests
# pass up to 1e-05 and for se(reading) in 2022 up to 1e-04
#
# In addition, we also display any individual discrepancies over 0.001 for
# diagnostic purposes, but these will not throw an error unless they cause the
# aggregate tolerance to be exceeded.
#
# If differences are observed only in a very small number of cases, the culprit
# is more likely to be small differences in the published PISA dataset and the
# dataset used to generate the official reports; these are usually identical but
# sometimes a handful of observations get added or disappear, presumably during
# last-minute quality control.
test_that("our estimates correspond to the official publication", {
  if (Sys.getenv("BRR_TEST_EQUIVALENCE") != "ON") skip()

  expected <- read_csv("../scores.csv", show_col_types = FALSE)
  expected$term <- str_c(expected$country_iso, ":", expected$cycle)

  pisa <- open_dataset(PISA_PATH) |>
    select(starts_with("pv") & ends_with("literacy"), "country_iso", "cycle", starts_with("w_literacy")) |>
    filter(cycle >= 2015, country_iso %in% OECD30_ISO) |>
    collect()
  fits <- brr(pv1literacy + pv2literacy + pv3literacy + pv4literacy + pv5literacy + pv6literacy + pv7literacy + pv8literacy + pv9literacy + pv10literacy ~ country_iso + cycle,
    statistic = weighted_mean_by,
    data = pisa,
    final_weights = matches("w_literacy_student_final"),
    replicate_weights = starts_with("w_literacy_student_r"),
    r = 80,
    .progress = FALSE
  )

  actual <- confint(fits, extra = TRUE)
  actual$domain <- "read"

  comparisons <- inner_join(actual, expected, by = c("term", "domain"))

  expect_equal(as.numeric(comparisons$estimate), comparisons$mean, tolerance = 1e-3)
  for (ix in which(abs(as.numeric(comparisons$estimate) - comparisons$mean) > 1e-3)) {
    warning(str_glue_data(comparisons[ix, ], "Discrepancy for {country_iso} {domain} in {cycle}: expected mu = {mean} but computed {estimate}"))
  }

  expect_equal(as.numeric(comparisons$se), comparisons$error, tolerance = 1e-3)
  for (ix in which(abs(as.numeric(comparisons$estimate) - comparisons$mean) > 1e-3)) {
    warning(str_glue_data(comparisons[ix, ], "Discrepancy for {country_iso} {domain} in {cycle}: expected se = {error} but computed {se}"))
  }
})

test_that("missing outcomes are handled appropriately", {
  pisa <- open_dataset(PISA_PATH) |>
    select(starts_with("pv") & ends_with("math"), "country_iso", "cycle", starts_with("w_math")) |>
    filter(cycle == 2000, country_iso == "BEL") |>
    collect()

  # there are two kinds of missingness here:
  # * we don't actually have pv6-10 (should be dealt with by `brr_bar`, `coef` and `confint`)
  # * math scores are not available for every student (should be dealt with by the `statistic`,
  #   but note that `brr` can pass an additional argument such as `na.rm` to `statistic`)
  fits10 <- brr(c("pv1math", "pv2math", "pv3math", "pv4math", "pv5math", "pv6math", "pv7math", "pv8math", "pv9math", "pv10math"),
    statistic = weighted_mean,
    data = pisa,
    final_weights = matches("w_math_student_final"),
    replicate_weights = starts_with("w_math_student_r"),
    .progress = FALSE,
    na_rm = TRUE
  )

  fits5 <- brr(c("pv1math", "pv2math", "pv3math", "pv4math", "pv5math"),
    statistic = weighted_mean,
    data = pisa,
    final_weights = matches("w_math_student_final"),
    replicate_weights = starts_with("w_math_student_r"),
    .progress = FALSE,
    na_rm = TRUE
  )

  # it is especially important to test `brr_var` because the `na_rm` logic
  # should make sure to count the non-NA outcomes instead of the outcomes from
  # the formula, otherwise the variance calculations will be off
  expect_true(!is.na(unname(coef(fits5))))
  expect_equal(unname(coef(fits10)), NaN)
  expect_equal(confint(fits10), tibble(term = "(Intercept)", estimate = NaN, lower = NaN, upper = NaN))
  expect_equal(coef(fits5), coef(fits10, na_rm = TRUE))
  expect_equal(brr_var(fits5), brr_var(fits10, na_rm = TRUE))
})
