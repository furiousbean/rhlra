prepare_cadzow_with_mgn <- function(obj, ...) UseMethod("prepare_cadzow_with_mgn")

prepare_cadzow_with_mgn.1d <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope) {
    if (is.null(series_for_cadzow)) series_for_cadzow <- series


    if (any(is.na(series))) {
        envelope <- correct_envelope(series, envelope)
    }

    left_diags <- inv_ac_diags(L, coefs)

    right_diags <- list(right_diag[effective_mask_for_weights(series, L)])

    list2env(prepare_cadzow(this, left_diags, right_diags), environment())

    common_expand_by_mask <- function(x) expand_by_mask(series, x)

    envelope_d <- Diagonal(length(series), 1/envelope)
    diags <- inv_ac_diags(length(series), coefs)

    ideal_weights <- envelope_d %*% band_mat_from_diags(diags) %*% envelope_d
    weights_chol <- chol(band_mat_from_diags(diags)) %*% envelope_d

    n <- length(series)
    df <- 2 * r + length(coefs) + 1

    list(series_for_cadzow = series_for_cadzow,
        envelope = envelope,
        left_diags = left_diags,
        left_mat = left_mat,
        left_chol = left_chol,
        left_chol_mat = left_chol_mat,
        right_diags = right_diags,
        right_mat = right_mat,
        right_chol = right_chol,
        right_chol_mat = right_chol_mat,
        weights = weights,
        common_expand_by_mask = common_expand_by_mask,
        ideal_weights = ideal_weights,
        weights_chol = weights_chol,
        n = n,
        df = df)
}

prepare_cadzow_with_mgn.1dm <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope) {
    if (is.null(series_for_cadzow)) series_for_cadzow <- series

    for (i in seq_along(series)) {
        if (any(is.na(series[[i]]))) {
            envelope[[i]] <- correct_envelope(series[[i]], envelope[[i]])
        }
    }

    left_diags <- sapply_ns(coefs, function(x) inv_ac_diags(L, x))

    right_diags <- sapply_ns(seq_along(series), function(i)
        list(right_diag[[i]][effective_mask_for_weights(series[[i]], L)]))

    list2env(prepare_cadzow(this, left_diags, right_diags), environment())

    common_expand_by_mask <- function(x) mapply(expand_by_mask, series, x, SIMPLIFY = FALSE)

    envelope_d <- sapply_ns(seq_along(series),
                         function(i) Diagonal(length(series[[i]]), 1/envelope[[i]]))
    diags <- sapply_ns(seq_along(series),
                    function(i) inv_ac_diags(length(series[[i]]), coefs[[i]]))
    ideal_weights <- sapply(seq_along(series), function(i)
        envelope_d[[i]] %*% band_mat_from_diags(diags[[i]]) %*% envelope_d[[i]])
    weights_chol <- sapply(seq_along(series), function(i)
        chol(band_mat_from_diags(diags[[i]])) %*% envelope_d[[i]])

    n <- sum(sapply(series, length))
    df <- (1 + length(series)) * r + length(coefs[[1]]) * (length(series) + 1)

    list(series_for_cadzow = series_for_cadzow,
        envelope = envelope,
        left_diags = left_diags,
        left_mat = left_mat,
        left_chol = left_chol,
        left_chol_mat = left_chol_mat,
        right_diags = right_diags,
        right_mat = right_mat,
        right_chol = right_chol,
        right_chol_mat = right_chol_mat,
        weights = weights,
        common_expand_by_mask = common_expand_by_mask,
        ideal_weights = ideal_weights,
        weights_chol = weights_chol,
        n = n,
        df = df)
}

cadzow_with_mgn <- function(this, series, L, r, coefs,
                            right_diag, series_for_cadzow = NULL,
                            use_mgn = TRUE, debug = FALSE,
                            envelope = unit_envelope(this, series), set_seed = NULL,
                            additional_pars, sylvester_nulling = NULL,
                            high_rank = FALSE) {
    series_for_mgn <- series

    list2env(prepare_cadzow_with_mgn(this, series, L, r, coefs,
                            right_diag, series_for_cadzow, envelope), environment())

    cadzow_call_list = list(this = this,
                            series = series_for_cadzow,
                            r = r,
                            left_chol_mat = left_chol_mat,
                            right_chol_mat = right_chol_mat,
                            weights_mat = weights,
                            set_seed = set_seed,
                            debug = debug,
                            sylvester_nulling = sylvester_nulling,
                            left_chol = left_chol,
                            right_chol = right_chol,
                            high_rank = high_rank)

    cadzow_data <- do.call(cadzow_iterations,
                           expand_pars_list(cadzow_call_list, cadzow_add_pars_names, additional_pars))

    answer <- cadzow_data

    answer$signal <- common_expand_by_mask(answer$signal)


    if (use_mgn) {
        if (high_rank) {
            stop("High rank with MGN is unsuppored")
        }

        mgn_call_list = list(this = this,
                             series = series_for_mgn,
                             signal = answer$signal,
                             r = r,
                             weights = ideal_weights,
                             debug = debug,
                             weights_chol = weights_chol)

        answer <- do.call(mgn, expand_pars_list(mgn_call_list, mgn_add_pars_names, additional_pars))

        answer$cadzow_it <- cadzow_data$it
    }

    if (!is.null(answer$noise_norm)) {
        answer$loglikelihood <- -n/2 * log(answer$noise_norm / n)
        answer$df <- df
    }

    answer
}

mgn_add_pars_names <- c("mgn_search_threshold", "mgn_it_limit", "mgn_j_limit", "mgn_max_step")
