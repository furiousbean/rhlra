hlra_igapfill <- function(series, r, igapfill_epsilon = 1e-6, igapfill_it_limit = 100, debug = TRUE,
                        set_seed = NULL, additional_pars = list()) {
    if (!is.null(set_seed)) {
        set_seed()
    }

    mask <- is.na(series)
    mask <- (1:length(series))[mask]

    signal <- series

    x <- rep(median(series, na.rm = TRUE), length(mask))
    prev <- NA
    series[mask] <- x
    dists <- numeric()
    dist <- 0
    step <- NA
    last_steps <- matrix(numeric(0), nrow = length(mask))

    it <- 0

    while (it < igapfill_it_limit && (anyNA(prev) || sum(abs(prev - x)) > igapfill_epsilon)) {
        pseudoobj <- list()
        class(pseudoobj) <- "1d"
        L <- default_L(series)

        series <- hlra_ssa(pseudoobj, series, L, r, debug, set_seed, additional_pars)$signal
        prev <- x
        x <- series[mask]
        series <- signal
        series[mask] <- x
        dist <- sum(abs(prev - x))
        dists <- c(dists, dist)
        it <- it + 1
    }

    if (debug == "P") {
        plot(series, type = "b", main = "igapfill: filled series", xlab = "index", ylab = "value")
        plot(dists, log = "y", type = "b", main = "igapfill: steps", xlab = "index", ylab = "value")
    }
    series
}

fill_gaps <- function(series, r, debug = FALSE, set_seed = NULL, additional_pars = list()) {
    result <- series
    if (any(is.na(series))) {
        if (debug == "P" || debug) cat("fillgaps: started\n")

        igapfill_call_list = list(series = series,
                                  r = r,
                                  set_seed = set_seed,
                                  debug = debug,
                                  additional_pars = additional_pars)

        result <- do.call(hlra_igapfill,
            expand_pars_list(igapfill_call_list, igapfill_add_pars_names, additional_pars))

        if (debug == "P" || debug) cat("fillgaps: done\n")
    }
    result
}

prepare_hlra_ssa <- function(obj, ...) UseMethod("prepare_hlra_ssa")

prepare_hlra_ssa.1d <- function(obj, series, L) {
    K <- length(series) - L + 1
    list(pseudord = rep(1, K), whitecoefs = numeric(0),
        pseudoenvelope = rep(1, length(series)))
}

prepare_hlra_ssa.1dm <- function(obj, series, L) {
    Ks <- sapply(series, function(series) length(as.numeric(series))) - L + 1
    list(pseudord = sapply(Ks, function(i) rep(1, i), simplify = FALSE),
        whitecoefs = sapply(seq_along(series), function(i) numeric(0), simplify = FALSE),
        pseudoenvelope = sapply(series, function(i) rep(1, i), simplify = FALSE))
}

hlra_ssa <- function(obj, series, L, r, debug, set_seed, additional_pars) {
    prepare_list <- prepare_hlra_ssa(obj, series, L)

    additional_pars["cadzow_it_limit"] <- 1
    additional_pars["cadzow_epsilon"] <- 1

    cadzow_with_mgn(obj, series, L, r, coefs = prepare_list$whitecoefs,
        right_diag = prepare_list$pseudord, use_mgn = FALSE, debug = debug,
        envelope = prepare_list$pseudoenvelope, series_for_cadzow = series,
        set_seed = set_seed, additional_pars = additional_pars)
}

igapfill_add_pars_names <- c("igapfill_epsilon", "igapfill_it_limit")
