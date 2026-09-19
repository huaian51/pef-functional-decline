# Shared helper functions for the public analysis package.

normalize_id <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", ".", "NA", "NaN", "NULL")] <- NA_character_
  sub("\\.0+$", "", z)
}

write_csv_utf8 <- function(x, path) {
  data.table::fwrite(x, path, bom = TRUE, na = "")
  message("Saved: ", path)
}

check_required <- function(dat, vars, label) {
  miss <- setdiff(vars, names(dat))
  if (length(miss)) {
    stop(label, " is missing required variable(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

usable_var <- function(x) {
  z <- x[!is.na(x)]
  length(z) > 0L && data.table::uniqueN(z) >= 2L
}

coverage_pct <- function(x) {
  100 * mean(!is.na(x))
}

cluster_vcov_glm <- function(fit, cluster) {
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  mu <- fitted(fit)

  pw <- weights(fit, type = "prior")
  if (is.null(pw)) pw <- rep(1, length(mu))

  W <- as.numeric(pw * mu * (1 - mu))
  bread0 <- crossprod(X, X * W)
  bread <- tryCatch(solve(bread0), error = function(e) MASS::ginv(bread0))

  score_i <- X * as.numeric(pw * (y - mu))
  score_g <- rowsum(score_i, group = as.factor(cluster), reorder = FALSE)
  meat <- crossprod(score_g)

  N <- nrow(X)
  p <- ncol(X)
  G <- nrow(score_g)
  correction <- if (G > 1L && N > p) {
    (G / (G - 1)) * ((N - 1) / (N - p))
  } else 1

  correction * bread %*% meat %*% bread
}

fit_cluster_logit <- function(dat, outcome, exposure, covars, cohort,
                              analysis_label, id_col = "id", weights_vec = NULL) {
  rhs <- c(exposure, covars)
  f <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  fit <- glm(f, data = as.data.frame(dat), family = binomial(), weights = weights_vec)
  V <- cluster_vcov_glm(fit, dat[[id_col]])

  b <- unname(coef(fit)[exposure])
  se <- sqrt(pmax(V[exposure, exposure], 0))

  data.table::data.table(
    cohort = cohort,
    analysis = analysis_label,
    interval_N = nrow(dat),
    person_N = data.table::uniqueN(dat[[id_col]]),
    event_N = sum(dat[[outcome]] == 1, na.rm = TRUE),
    beta = b,
    SE = se,
    OR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se),
    P = 2 * pnorm(abs(b / se), lower.tail = FALSE)
  )
}

build_estimand_data <- function(dat, estimand) {
  d <- data.table::copy(dat)
  d[, destination := as.character(destination)]

  if (estimand == "A_ANY_LIMITATION") {
    d <- d[destination %in% c("Independent", "Mild limitation", "Severe limitation")]
    d[, y := as.integer(destination != "Independent")]
    return(d)
  }
  if (estimand == "MILD_VS_INDEPENDENT") {
    d <- d[destination %in% c("Independent", "Mild limitation")]
    d[, y := as.integer(destination == "Mild limitation")]
    return(d)
  }
  if (estimand == "SEVERE_VS_INDEPENDENT") {
    d <- d[destination %in% c("Independent", "Severe limitation")]
    d[, y := as.integer(destination == "Severe limitation")]
    return(d)
  }
  if (estimand == "SEVERE_VS_MILD") {
    d <- d[destination %in% c("Mild limitation", "Severe limitation")]
    d[, y := as.integer(destination == "Severe limitation")]
    return(d)
  }
  stop("Unknown estimand: ", estimand, call. = FALSE)
}

fit_multinom <- function(dat, exposure = "burden_z", covars) {
  f <- as.formula(paste("dest3 ~", paste(c(exposure, covars), collapse = " + ")))
  nnet::multinom(
    f,
    data = as.data.frame(dat),
    trace = FALSE,
    maxit = 1000,
    MaxNWts = 20000,
    reltol = 1e-10
  )
}

extract_multinom_exposure <- function(fit, exposure = "burden_z") {
  cf <- coef(fit)
  if (is.vector(cf)) stop("Multinomial coefficient matrix collapsed unexpectedly.")
  need_rows <- c("Mild limitation", "Severe limitation")
  if (!all(need_rows %in% rownames(cf))) {
    stop("Unexpected multinomial outcome rows: ", paste(rownames(cf), collapse = ", "))
  }
  if (!(exposure %in% colnames(cf))) stop("Exposure coefficient not found in multinomial model.")

  b_mild <- unname(cf["Mild limitation", exposure])
  b_severe <- unname(cf["Severe limitation", exposure])
  c(mild = b_mild, severe = b_severe, contrast = b_severe - b_mild)
}

bootstrap_multinom_person <- function(dat, exposure = "burden_z", covars,
                                      B = 1000L, seed = 1L, progress_every = 100L) {
  set.seed(seed)
  ids_chr <- as.character(unique(dat$id))
  n_ids <- length(ids_chr)
  if (n_ids < 100L) stop("Too few participants for cluster bootstrap.")

  out <- matrix(
    NA_real_, nrow = B, ncol = 3,
    dimnames = list(NULL, c("beta_mild", "beta_severe", "beta_severe_minus_mild"))
  )

  f <- as.formula(paste("dest3 ~", paste(c(exposure, covars), collapse = " + ")))
  id_row_chr <- as.character(dat$id)

  for (b in seq_len(B)) {
    sampled <- sample(ids_chr, size = n_ids, replace = TRUE)
    mult <- table(sampled)
    keep_ids <- names(mult)
    keep <- id_row_chr %in% keep_ids
    db <- data.table::copy(dat[keep])
    db[, boot_w := as.numeric(mult[as.character(id)])]

    fit <- tryCatch(
      nnet::multinom(
        f, data = as.data.frame(db), weights = db$boot_w,
        trace = FALSE, maxit = 1000, MaxNWts = 20000, reltol = 1e-10
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      z <- tryCatch(extract_multinom_exposure(fit, exposure), error = function(e) NULL)
      if (!is.null(z) && all(is.finite(z))) out[b, ] <- z
    }

    if (progress_every > 0L && b %% progress_every == 0L) {
      message("  bootstrap ", b, "/", B)
    }
  }
  out
}

meta_reml_mkh <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  yi <- beta[ok]
  sei <- se[ok]
  vi <- sei^2
  k <- length(yi)

  if (k < 2L) {
    return(data.table::data.table(
      k = k, pooled_beta = NA_real_, pooled_OR = NA_real_, pooled_effect = NA_real_,
      CI_low_mKH = NA_real_, CI_high_mKH = NA_real_, P_mKH = NA_real_,
      tau2_REML = NA_real_, I2 = NA_real_, Q = NA_real_, q_HK = NA_real_,
      q_modified = NA_real_
    ))
  }

  wf <- 1 / vi
  mu_fixed <- sum(wf * yi) / sum(wf)
  Q <- sum(wf * (yi - mu_fixed)^2)
  I2 <- if (Q > 0) max(0, 100 * (Q - (k - 1L)) / Q) else 0

  reml_nll <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }

  vyi <- suppressWarnings(var(yi))
  upper <- max(1, ifelse(is.finite(vyi), vyi * 100, 1))
  tau2 <- optimize(reml_nll, interval = c(0, upper))$minimum
  wr <- 1 / (vi + tau2)
  mu <- sum(wr * yi) / sum(wr)
  q_hk <- sum(wr * (yi - mu)^2) / (k - 1L)
  q_mkh <- max(1, q_hk)
  se_mkh <- sqrt(q_mkh / sum(wr))
  tcrit <- qt(0.975, df = k - 1L)
  p <- 2 * pt(abs(mu / se_mkh), df = k - 1L, lower.tail = FALSE)

  data.table::data.table(
    k = k,
    pooled_beta = mu,
    pooled_OR = exp(mu),
    pooled_effect = exp(mu),
    CI_low_mKH = exp(mu - tcrit * se_mkh),
    CI_high_mKH = exp(mu + tcrit * se_mkh),
    P_mKH = p,
    tau2_REML = tau2,
    I2 = I2,
    Q = Q,
    q_HK = q_hk,
    q_modified = q_mkh
  )
}

run_meta_by_analysis <- function(results, analysis_col = "analysis") {
  groups <- unique(results[[analysis_col]])
  out <- lapply(groups, function(g) {
    z <- results[get(analysis_col) == g]
    m <- meta_reml_mkh(z$beta, z$SE)
    m[, analysis := g]
    m[, positive_cohort_N := sum(z$OR > 1, na.rm = TRUE)]
    m
  })
  data.table::rbindlist(out, fill = TRUE)
}
