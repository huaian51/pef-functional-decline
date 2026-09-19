# Synthetic demonstration data ONLY.
# These data are not derived from HRS, ELSA, SHARE, or CHARLS and must not be
# used for scientific interpretation. They exist only to demonstrate the input
# schema and allow a smoke test of the public code.

rm(list = ls())
set.seed(20260919)
suppressPackageStartupMessages(library(data.table))

dir.create("data_harmonized", recursive = TRUE, showWarnings = FALSE)
cohorts <- c("HRS", "ELSA", "SHARE", "CHARLS")

for (j in seq_along(cohorts)) {
  coh <- cohorts[j]
  N <- 420L + 20L * j
  id <- sprintf("%s_%05d", coh, seq_len(N))
  sex <- sample(c("Male", "Female"), N, replace = TRUE)
  height <- rnorm(N, ifelse(sex == "Male", 172, 160), 7)
  age1 <- rnorm(N, 66 + j, 7)
  age2 <- age1 + 4
  pef1 <- 520 - 3.0 * age1 + 1.7 * height + ifelse(sex == "Male", 55, 0) + rnorm(N, 0, 45)
  latent <- rnorm(N)
  pef2 <- 520 - 3.0 * age2 + 1.7 * height + ifelse(sex == "Male", 55, 0) + 10 * latent + rnorm(N, 0, 45)

  pef <- data.table(
    id = id,
    pef_1 = pef1, age_1 = age1, sex_1 = sex, height_1 = height,
    pef_2 = pef2, age_2 = age2, sex_2 = sex, height_2 = height
  )
  saveRDS(pef, file.path("data_harmonized", paste0(coh, "_pef_exposure.rds")))

  # Two possible functional intervals per participant.
  intervals <- rbindlist(lapply(1:2, function(k) {
    age <- age2 + 2 * (k - 1)
    bmi <- rnorm(N, 26, 4)
    edu <- sample(0:3, N, replace = TRUE)
    smoke <- rbinom(N, 1, 0.18)
    risk <- plogis(-2.5 + 0.25 * latent + 0.02 * (age - 65) + 0.15 * smoke)
    any <- rbinom(N, 1, risk)
    severe <- ifelse(any == 1, rbinom(N, 1, plogis(-0.7 + 0.35 * latent)), 0)
    dest <- ifelse(any == 0, "Independent", ifelse(severe == 1, "Severe limitation", "Mild limitation"))
    obs <- rbinom(N, 1, plogis(2.8 - 0.15 * latent))
    dest[obs == 0] <- NA_character_

    data.table(
      id = id,
      interval_f = paste0("interval_", k),
      destination = dest,
      age = age,
      sex = sex,
      height = height,
      bmi = bmi,
      education = factor(edu),
      observed_next = obs,
      smoking = smoke,
      ever_smoking = pmax(smoke, rbinom(N, 1, 0.25)),
      lung_disease = rbinom(N, 1, 0.08),
      diabetes = rbinom(N, 1, 0.16),
      hypertension = rbinom(N, 1, 0.42),
      heart_disease = rbinom(N, 1, 0.12),
      stroke = rbinom(N, 1, 0.06),
      cancer = rbinom(N, 1, 0.08),
      depression = rnorm(N),
      adl_iadl_t0 = 0L,
      adl_iadl_t1 = ifelse(is.na(dest), NA_integer_, ifelse(dest == "Independent", 0L, ifelse(dest == "Mild limitation", 1L, sample(2:4, N, replace = TRUE))))
    )
  }))
  saveRDS(intervals, file.path("data_harmonized", paste0(coh, "_functional_intervals.rds")))
}

cat("Synthetic demonstration inputs written to ./data_harmonized/.\n")
