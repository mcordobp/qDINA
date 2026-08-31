############################################################
# qDINA empirical application
############################################################
# libraries


library(Rcpp)
library(RcppArmadillo)
library(rgen)
library(dina)
library(simcdm)
library(CDM)
library(dplyr)
library(tidyr)


############################################################
# Load qDINA implementation and phi calibration
############################################################
Rcpp::sourceCpp("../scr/qdina_gibbs.cpp")
source("../R/phi_calibration.R")

############################################################
# Load empirical data
############################################################

data(data.cdm05, package = "CDM")
dat <- data.cdm05$data
Q <- data.cdm05$q.matrix

############################################################
# Fit DINA and qDINA
############################################################
system.time({
  # Baseline DINA model
  run_DINA_base <- dina(
    Y = as.matrix(dat),
    Q = as.matrix(Q),
    chain_length = 20000
  )

  # estimate item-specific dispersion parameters
  out_phi <- phi_BBQ_window(
    run_DINA_base,
    as.matrix(dat),
    as.matrix(Q),
    t_start = 2000,
    t_end = 20000,
    weight = "uniform"
  )

  phi_bbq_final <- pmax(1, out_phi$phi_final)
  phi_bbq_mean  <- pmax(1, out_phi$phi_mean)

  # qDINA using thresholded phi estimates
  run_qDINA_Gibbs_final <- DINA_Gibbs_cpp_phi_item(
    Y = as.matrix(dat),
    Q = as.matrix(Q),
    phi = phi_bbq_final,
    chain_length = 20000
  )

  # qDINA using posterior mean phi estimates
  run_qDINA_Gibbs_mean <- DINA_Gibbs_cpp_phi_item(
    Y = as.matrix(dat),
    Q = as.matrix(Q),
    phi = phi_bbq_mean,
    chain_length = 20000
  )
})

############################################################
# item-specific dispersion summary
############################################################

table_phi <- data.frame(
  item = seq_along(out_phi$phi_mean),
  phi_mean = round(out_phi$phi_mean, 3),
  phi_median = round(out_phi$phi_median, 3),
  p_phi_gt_1 = round(out_phi$p_phi_gt, 3),
  phi_p025 = round(
    apply(
      out_phi$phi_by_iter,
      1,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    3
  ),
  phi_p975 = round(
    apply(
      out_phi$phi_by_iter,
      1,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    ),
    3
  )
)

table_phi

############################################################
# latent class proportion summaries
############################################################

summarize_PIs <- function(fit, model_name,
                          t_start = 2000, t_end = 20000, by = 2) {

  iter_idx <- seq(t_start, t_end, by = by)
  pis_sub <- fit$PIs[, iter_idx, drop = FALSE]

  data.frame(
    model = model_name,
    class = 0:(nrow(pis_sub) - 1),
    pi_mean = round(100 * apply(pis_sub, 1, mean, na.rm = TRUE), 1),
    pi_p025 = round(100 * apply(pis_sub, 1, quantile,
                                probs = 0.025, na.rm = TRUE), 1),
    pi_p975 = round(100 * apply(pis_sub, 1, quantile,
                                probs = 0.975, na.rm = TRUE), 1)
  )
}

table_pi <- rbind(
  summarize_PIs(run_DINA_base, "DINA"),
  summarize_PIs(run_qDINA_Gibbs_mean, "qDINA")
)

table_pi

############################################################
# wide latent class table
############################################################

table_pi_wide <- table_pi %>%
  pivot_wider(
    names_from = model,
    values_from = c(
      pi_mean,
      pi_p025,
      pi_p975
    )
  )

table_pi_wide
