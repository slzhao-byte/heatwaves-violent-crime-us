# Heatwaves and police-recorded violent offenses in US cities: association,
# burden attributable to anthropogenic warming and projected burden.
#
# Runs on example_data.rds, simulated data for four hypothetical cities.
# Study inputs with the same tables and columns can replace it (README).
# Results are written to output/.

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
dir.create("output", showWarnings = FALSE)

inputs <- lapply(readRDS("example_data.rds"), as.data.table)
daily <- inputs$daily           # city-day weather, holidays and offenses
cities <- inputs$cities         # overall Social Vulnerability Index
population <- inputs$population # official city population, 2005-2022

years <- 2005:2022
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
groups <- c("low", "high")

# Heatwave day: a run of at least two days at or above the threshold,
# identified on a complete daily calendar.
heatwave <- function(temperature, threshold) {
  runs <- rle(temperature >= threshold)
  rep(as.integer(runs$values & runs$lengths >= 2L), runs$lengths)
}

# Heatwave history of days t-6 to t coded as 0-127 (lag 0 in bit 0).
lag_pattern <- function(flags) {
  as.integer(do.call(cbind, shift(flags, 0:6)) %*% 2^(0:6))
}

# Sum of weights by year (rows) and lag history (columns).
pattern_weights <- function(year, pattern, weight, all_years) {
  sums <- data.table(year, pattern, weight)[
    , list(weight = sum(weight)), by = list(year, pattern)
  ]
  out <- matrix(0, length(all_years), 128L)
  out[cbind(match(sums$year, all_years), sums$pattern + 1L)] <- sums$weight
  out
}

# Estimate (first element) and 95% CI from the joint draws (the rest).
ci <- function(x) {
  x <- unname(x)
  draws <- x[-1L]
  q <- quantile(draws[is.finite(draws)], c(0.025, 0.975), names = FALSE)
  c(estimate = x[1L], lower = q[1L], upper = q[2L])
}

# 1. Association ------------------------------------------------------------
# Time-stratified case-crossover design, conditional Poisson model with a
# distributed lag term for heatwave days over lags 0-6 and lag 0-7 terms for
# relative humidity and precipitation.
setorder(daily, city, date)
daily[, `:=`(year = year(date), month = month(date), dow = wday(date))]
daily[, warm := month %in% 5:9 & year %in% years]
cities[, group := ifelse(svi >= quantile(svi, 0.75), "high", "low")]

lag_args <- list(fun = "ns", knots = dlnm::logknots(6, 2), intercept = TRUE)
lag_basis <- do.call(dlnm::onebasis, c(list(x = 0:6), lag_args))
rh_var <- dlnm::onebasis(daily$humidity, fun = "ns", df = 3)
rain_var <- dlnm::onebasis(daily$precipitation, fun = "ns", df = 3)
rh_arg <- list(fun = "ns", knots = attr(rh_var, "knots"),
               Boundary.knots = attr(rh_var, "Boundary.knots"))
rain_arg <- list(fun = "ns", knots = attr(rain_var, "knots"),
                 Boundary.knots = attr(rain_var, "Boundary.knots"))
hw_terms <- paste0("hw", 1:4)
rh_terms <- paste0("rh", 1:9)
rain_terms <- paste0("rain", 1:9)

panel <- list()
for (id in cities$city) {
  d <- daily[city == id]
  hw <- heatwave(d$temperature, quantile(d$temperature[d$warm], 0.9))
  cb <- dlnm::crossbasis(hw, lag = 6, argvar = list(fun = "lin"),
                         arglag = lag_args)
  rh <- dlnm::crossbasis(d$humidity, lag = 7, argvar = rh_arg,
                         arglag = list(fun = "ns", df = 3))
  rain <- dlnm::crossbasis(d$precipitation, lag = 7, argvar = rain_arg,
                           arglag = list(fun = "ns", df = 3))
  d[, c(hw_terms, rh_terms, rain_terms) := as.data.table(cbind(cb, rh, rain))]
  panel[[id]] <- d[warm & !is.na(offenses)]
}
panel <- rbindlist(panel)
panel[, stratum := interaction(city, year, month, dow, drop = TRUE)]
panel <- panel[complete.cases(panel[, c(hw_terms, rh_terms, rain_terms),
                                    with = FALSE])]
panel <- panel[, if (sum(offenses) > 0) .SD, by = stratum]
high <- cities$group[match(panel$city, cities$city)] == "high"
panel[, paste0("hw_high", 1:4) := lapply(.SD, function(x) x * high),
      .SDcols = hw_terms]

covariates <- paste(c("holiday", rh_terms, rain_terms), collapse = " + ")
model_terms <- list(all_cities = hw_terms,
                    by_group = c(hw_terms, paste0("hw_high", 1:4)))
models <- list()
for (m in names(model_terms)) {
  terms <- model_terms[[m]]
  fit <- fixest::fepois(
    as.formula(paste("offenses ~", paste(terms, collapse = " + "), "+",
                     covariates, "| stratum")),
    data = panel, notes = FALSE
  )
  v <- vcov(fit, vcov = "iid",
            ssc = fixest::ssc(adj = FALSE, fixef.K = "none"))
  models[[m]] <- list(b = coef(fit)[terms], V = v[terms, terms])
}

# Lag-specific and cumulative (lags 0-6) relative risks.
z <- qnorm(0.975)
basis <- rbind(lag_basis, colSums(lag_basis))
zero <- matrix(0, nrow(basis), 4L)
contrasts <- list(all = list("all_cities", basis),
                  low = list("by_group", cbind(basis, zero)),
                  high = list("by_group", cbind(basis, basis)))
association <- list()
for (g in names(contrasts)) {
  m <- models[[contrasts[[g]][[1L]]]]
  k <- contrasts[[g]][[2L]]
  log_rr <- drop(k %*% m$b)
  se <- sqrt(rowSums((k %*% m$V) * k))
  association[[g]] <- data.table(
    group = g, lag = c(0:6, "0-6"), rr = exp(log_rr),
    lower = exp(log_rr - z * se), upper = exp(log_rr + z * se)
  )
}
# Ratio of cumulative RRs (high vs low) and two-sided Wald test.
k <- c(rep(0, 4L), colSums(lag_basis))
log_ratio <- sum(k * models$by_group$b)
se <- sqrt(drop(k %*% models$by_group$V %*% k))
association$ratio <- data.table(
  group = "high_vs_low", lag = "0-6", rr = exp(log_ratio),
  lower = exp(log_ratio - z * se), upper = exp(log_ratio + z * se),
  p_interaction = pchisq((log_ratio / se)^2, 1L, lower.tail = FALSE)
)
fwrite(rbindlist(association, fill = TRUE), "output/association.csv")

# 2. Local warming and counterfactual temperatures --------------------------
# Smoothed HadCRUT5 GMST anomalies (central, lower and upper series), city
# warming per unit GMST from warm-season ERA5-Land means for 1950-2022, and
# nine counterfactual daily series (three GMST series x slope and its 95% CI).
gmst <- inputs$gmst
setorder(gmst, year, month)
local <- inputs$local_monthly
stopifnot(local[, .N, by = .(city, year)][, all(N == 5L)])
local_annual <- local[, .(temperature = mean(temperature)), by = .(city, year)]
baseline <- gmst$year %in% 1850:1900
variants <- character()
slopes <- list()
for (s in c("central", "lower", "upper")) {
  x <- gmst[[s]] - mean(gmst[[s]][baseline])
  ssa <- Rssa::ssa(x, L = 120L, neig = 120L, kind = "1d-ssa",
                   svd.method = "eigen")
  trend <- as.numeric(Rssa::reconstruct(ssa, groups = list(1L))[[1L]])
  gmst[, (s) := trend - mean(trend[baseline])]
  annual <- gmst[year %in% 1950:2022, .(anomaly = mean(get(s))), by = year]
  slopes[[s]] <- merge(local_annual, annual, by = "year")[, {
    fit <- lm(temperature ~ anomaly)
    bounds <- confint(fit, "anomaly")
    .(fit = coef(fit)[["anomaly"]], lower = bounds[1L], upper = bounds[2L])
  }, by = city]
  if (s == "central") {
    warming <- CJ(city = slopes$central$city, year = annual$year)
    warming[, warming := slopes$central$fit[match(city, slopes$central$city)] *
              annual$anomaly[match(year, annual$year)]]
    fwrite(warming, "output/local_warming.csv")
  }
  anomaly <- gmst[daily, on = .(year, month), get(s)]
  for (b in c("fit", "lower", "upper")) {
    kappa <- slopes[[s]][[b]][match(daily$city, slopes[[s]]$city)]
    v <- paste(s, b, sep = "_")
    daily[, (v) := temperature - ifelse(warm, kappa * anomaly, 0)]
    variants <- c(variants, v)
  }
}

# 3. Burden attributable to anthropogenic warming ---------------------------
# Backward attributable numbers in the factual climate (observed temperature)
# and the counterfactual climates, with each group's exposure-response
# function; 10,000 joint coefficient draws and the nine climate variants.
worlds <- c("temperature", variants)
histories <- sapply(0:6, function(l) bitwAnd(bitwShiftR(0:127, l), 1L))
design <- histories %*% lag_basis
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261110L)
b <- models$by_group$b
draws <- rbind(b, MASS::mvrnorm(10000L, b, models$by_group$V))
theta <- list(low = draws[, 1:4], high = draws[, 1:4] + draws[, 5:8])
af <- lapply(theta, function(x) -expm1(-design %*% t(x)))

summarize_burden <- function(totals, person_years, recorded) {
  factual <- totals[1L, ]
  delta <- sweep(-totals[-1L, , drop = FALSE], 2L, factual, "+")
  proportion <- sweep(delta, 2L, factual, "/")
  proportion[, factual <= 0 & seq_along(factual) > 1L] <- NA_real_
  pool <- function(x) c(x[1L, 1L], as.vector(x[, -1L, drop = FALSE]))
  out <- rbind(
    factual_offenses = ci(factual),
    counterfactual_offenses = ci(pool(totals[-1L, , drop = FALSE])),
    attributable_offenses = ci(pool(delta)),
    attributable_proportion = ci(pool(proportion)),
    per_million_person_years = ci(1e6 * pool(delta) / person_years),
    share_of_recorded_offenses = ci(pool(delta) / recorded)
  )
  as.data.table(out, keep.rownames = "metric")
}

crime_years <- daily[warm & !is.na(offenses),
                     .(offenses = sum(offenses)), by = .(city, year)]
support <- merge(crime_years, population, by = c("city", "year"))
weights <- lapply(setNames(nm = groups), function(g) {
  array(0, c(length(years), 128L, length(worlds)))
})
person_years <- recorded <- lapply(setNames(nm = groups), function(g) {
  numeric(length(years))
})
attribution <- list()
for (id in cities$city) {
  d <- daily[city == id]
  g <- cities[city == id, group]
  take <- d$warm & !is.na(d$offenses)
  threshold <- quantile(d$temperature[d$warm], 0.9)
  w <- array(0, c(length(years), 128L, length(worlds)))
  for (k in seq_along(worlds)) {
    hw <- heatwave(d[[worlds[k]]], threshold)
    w[, , k] <- pattern_weights(d$year[take], lag_pattern(hw)[take],
                                d$offenses[take], years)
  }
  s <- support[city == id]
  j <- match(s$year, years)
  weights[[g]] <- weights[[g]] + w
  person_years[[g]][j] <- person_years[[g]][j] + s$population
  recorded[[g]][j] <- recorded[[g]][j] + s$offenses
  totals <- t(colSums(w)) %*% af[[g]]
  attribution[[id]] <- summarize_burden(
    totals, sum(s$population), sum(s$offenses)
  )[, `:=`(level = id, year = "2005-2022")]
}

period <- factual <- central <- list()
for (g in groups) {
  period[[g]] <- t(colSums(weights[[g]])) %*% af[[g]]
  factual[[g]] <- weights[[g]][, , 1L] %*% af[[g]]
  central[[g]] <- weights[[g]][, , 2L] %*% af[[g]]
}
period$all <- period$low + period$high
factual$all <- factual$low + factual$high
central$all <- central$low + central$high
person_years$all <- person_years$low + person_years$high
recorded$all <- recorded$low + recorded$high
for (g in c(groups, "all")) {
  attribution[[g]] <- summarize_burden(
    period[[g]], sum(person_years[[g]]), sum(recorded[[g]])
  )[, `:=`(level = g, year = "2005-2022")]
  # Annual estimates use the central GMST series and fitted city slopes.
  for (i in seq_along(years)) {
    attribution[[paste(g, years[i])]] <- summarize_burden(
      rbind(factual[[g]][i, ], central[[g]][i, ]),
      person_years[[g]][i], recorded[[g]][i]
    )[, `:=`(level = g, year = as.character(years[i]))]
  }
}
delta <- lapply(period, function(p) {
  sweep(-p[-1L, , drop = FALSE], 2L, p[1L, ], "+")
})
for (g in groups) {
  share <- delta[[g]] / delta$all
  share[delta$all <= 0 & col(share) > 1L] <- NA_real_
  attribution[[paste(g, "share")]] <- data.table(
    metric = "share_of_attributable_offenses",
    t(ci(c(share[1L, 1L], share[, -1L]))), level = g, year = "2005-2022"
  )
}
attribution <- rbindlist(attribution, use.names = TRUE)
fwrite(attribution, "output/attribution.csv")

# Linear trend in the annual proportion attributable (points per decade).
annual <- attribution[metric == "attributable_proportion" &
                        year != "2005-2022"]
trend <- annual[, {
  fit <- lm(1000 * estimate ~ as.numeric(year))
  bounds <- confint(fit)[2L, ]
  .(points_per_decade = coef(fit)[[2L]], lower = bounds[[1L]],
    upper = bounds[[2L]])
}, by = level]
fwrite(trend, "output/attribution_trend.csv")

# 4. Projected heatwave-related violent offenses ----------------------------
# Equation (3): non-heatwave offense rate x population x (RR - 1) for the
# projected heatwave history; thresholds fixed at the 2005-2022 90th
# percentile of each scenario's bias-corrected temperatures; 1,000 draws.
projected <- inputs$projected_temperature
projected_population <- inputs$projected_population
setkey(projected, city, scenario, date)
setkey(projected_population, city, scenario, year)
all_years <- 2005:2100
calendar <- seq(as.Date("2005-01-01"), as.Date("2100-12-31"), by = "day")
calendar_year <- year(calendar)
calendar_warm <- month(calendar) %in% 5:9
reference <- calendar_warm & calendar_year %in% years

set.seed(20261110L)
draws <- rbind(b, MASS::mvrnorm(1000L, b, models$by_group$V))
theta <- list(low = draws[, 1:4], high = draws[, 1:4] + draws[, 5:8])
rr_minus_1 <- lapply(theta, function(x) expm1(design %*% t(x)))

expected <- people <- list()
for (s in scenarios) for (g in groups) {
  expected[[s]][[g]] <- matrix(0, length(all_years), 128L)
  people[[s]][[g]] <- numeric(length(all_years))
}
for (id in cities$city) {
  x <- daily[city == id]
  g <- cities[city == id, group]
  hw <- heatwave(x$temperature, quantile(x$temperature[x$warm], 0.9))
  quiet <- x[warm & !is.na(offenses) & hw == 0L]
  official <- population[city == id]
  rate <- sum(quiet$offenses) /
    sum(official$population[match(quiet$year, official$year)])
  for (s in scenarios) {
    d <- projected[.(id, s)]
    pop <- projected_population[.(id, s)]$population
    stopifnot(identical(d$date, calendar), length(pop) == length(all_years))
    hw <- heatwave(d$temperature, quantile(d$temperature[reference], 0.9))
    w <- pattern_weights(calendar_year[calendar_warm],
                         lag_pattern(hw)[calendar_warm],
                         rep(1, sum(calendar_warm)), all_years)
    expected[[s]][[g]] <- expected[[s]][[g]] + rate * pop * w
    people[[s]][[g]] <- people[[s]][[g]] + pop
  }
}

windows <- data.table(
  window = c("2005-2022", "2050-2059", "2090-2099",
             paste0("moving10_", 2005:2091)),
  start = c(2005L, 2050L, 2090L, 2005:2091),
  end = c(2022L, 2059L, 2099L, 2014:2100)
)
ref <- match(years, all_years)
projection <- list()
for (s in scenarios) {
  number <- lapply(groups, function(g) expected[[s]][[g]] %*% rr_minus_1[[g]])
  names(number) <- groups
  number$all <- number$low + number$high
  pyears <- people[[s]]
  pyears$all <- pyears$low + pyears$high
  rate_at <- function(g, j) {
    1e6 * colSums(number[[g]][j, , drop = FALSE]) / sum(pyears[[g]][j])
  }
  for (i in seq_len(nrow(windows))) {
    j <- match(windows$start[i]:windows$end[i], all_years)
    gap <- rate_at("high", j) - rate_at("low", j)
    rows <- list(
      rate_difference_high_minus_low = gap,
      difference_change_from_reference = gap -
        (rate_at("high", ref) - rate_at("low", ref))
    )
    for (g in c(groups, "all")) {
      n <- colSums(number[[g]][j, , drop = FALSE])
      rows[[paste0(g, ":offenses")]] <- n
      rows[[paste0(g, ":annual_offenses")]] <- n / length(j)
      rows[[paste0(g, ":rate_per_million_person_years")]] <- rate_at(g, j)
      rows[[paste0(g, ":rate_change_from_reference")]] <-
        rate_at(g, j) - rate_at(g, ref)
    }
    projection[[length(projection) + 1L]] <- data.table(
      scenario = s, window = windows$window[i], metric = names(rows),
      do.call(rbind, lapply(rows, ci))
    )
  }
}
fwrite(rbindlist(projection), "output/projection.csv")
cat("Done: results in output/.\n")
