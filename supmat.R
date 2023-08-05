# Implementation of Ser-Giacomi et al. fit (ref. 29)
s <- function(a, b, r, s, ploc) {
  sloc <-  tryCatch({
    h1 <- BAS::hypergeometric2F1(1, 1 + a, 2 + b,
                                 exp(-r) * (1 - ploc),
                                 log = FALSE)
    h2 <- BAS::hypergeometric2F1(1, 1 + a, 2 + b, exp(-r), log = FALSE)
    res <- s * (1 - ((1 - ploc) * h1 / h2))
    msg <- "no issues"
    list(predict = res, predict.debug = msg)
  },
  warning = function(cond) {
    res <- NA
    msg <-  as.character(cond)
    list(predict = res, predict.debug = msg)
  }
  )
  return(sloc)
}

#needed by batch_downscale, calculates adjusted total abundance and richness
#after removing all the species above threshold.
imp2 <- function(dat, thresh) {
  dat %>%
    group_by(station) %>%
    filter(Count < thresh) %>%
    summarise(tot = sum(Count), rich = n(), .groups = "drop")
}

#needed by batch_downscale, calculates the p needed by the estimator
imp3 <- function(refst, counts) {
  refcount <- filter(counts, station == refst) %>%
    pull(tot)
  res <- dplyr::mutate(counts, p = tot / refcount)
  return(res)
}

# implementation of the batch downscaling function detailed in the paper
batch_downscale <- function(x) {
  fits <- x %>%
    rename(station.ref = station) %>%
    mutate(adj.tc.r = map(vec_x_max, ~imp2(ext_data, .x))) %>%
    nest(fit_parms = -c(station.ref, adj.tc.r)) %>%
    mutate(p = map2(station.ref, adj.tc.r, ~imp3(.x, .y))) %>%
    unnest(p) %>%
    filter(p < 1) %>%
    mutate(predict = map2(fit_parms, p, ~s(.x$vec_estimated_l,
                                           .x$vec_estimated_k,
                                           .x$vec_estimated_r,
                                           .x$vec_n, .y)
    )) %>%
    unnest_wider(predict)
  return(fits)
}

library(tidyverse)
load("all_best.Rdata") #fit parameters for all stations (fits with no xmax_tsh)
load("ext_data.Rdata")

#Calculate downscaling using each station as reference, for all stations whose richness is lower than the reference
batch_ds <- batch_downscale(all_best)

# remove all inconclusive stations
wt <- batch_ds %>%
  unnest(fit_parms) %>%
  filter(vec_n > rich, rich > 10) %>%
  select(-adj.tc.r)


# create an expanded and annotate object for plotting
wt_w_spe <- wt %>%
  mutate(spe = (predict - rich) / rich) %>%
  nest(data = -station.ref) %>%
  rename(station = station.ref) %>%
  rename(station.ref = station) %>%
  unnest(data) %>%
  select(station.ref, station, rich, predict, spe) %>%
  separate(station, into = c("stat", "Depth"), remove = FALSE) %>%
  add_count(station.ref) %>%
  filter(grepl("SRF", station.ref), Depth == "SRF") %>%
  mutate(rich_cl = cut(rich, c(10, 18, 32, 58, 104, 186, 334, 600)))

#base plot
plo <- ggplot(wt_w_spe) +
  geom_vline(xintercept = "173_SRF", lty = 2, col = "red", lwd = 0.2) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = c(-1, 12)) +
  scale_fill_steps("% of predicted stations", low = "black",
                   high = "dodgerblue") +
  ggthemes::theme_tufte(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))   +
  labs(caption = "64 (n=291) and 135 (n=139)
                  are cropped out (-5.92 and -2.19 average)")

#boxplot logplot
plo +
  aes(x = rich_cl, y = spe) +
  xlab("Reference station richness class") +
  ylab("\U0394 S")
