############################################################
# qDINA simulation fitting script
############################################################
# Libraries

library(Rcpp)
library(RcppArmadillo)
library(rgen)
library(dina)
library(simcdm)
library(parallel)

############################################################
# compile qDINA Gibbs sampler
############################################################

cpp_file <- "./src/qdina_gibbs.cpp"
Rcpp::sourceCpp(cpp_file)


############################################################
# item-specific dispersion calibration
############################################################

phi_BBQ_window <- function(fit, Y, Q,
                           t_start = 2000, t_end = 20000,
                           weight = c("uniform", "pi", "nc_over_N"),
                           eps = 1e-12,
                           prob_cut = 0.95,
                           delta = 0,
                           use_stat = c("mean", "median")) {

  weight   <- match.arg(weight)
  use_stat <- match.arg(use_stat)

  CL       <- fit$CLASSES
  ss_chain <- fit$SigS
  gs_chain <- fit$GamS
  pi_chain <- fit$PIs

  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  C <- 2^K
  # Determine whether latent classes are indexed from 0 or 1
  cls_min <- min(CL[, t_start:t_end], na.rm = TRUE)
  offset  <- if (cls_min == 0) 1 else 0
  # Attribute patterns and ideal-response indicators
  Amat <- simcdm::attribute_classes(K)
  AQ   <- Amat %*% t(Q)
  qq   <- rowSums(Q * Q)

  Tidx    <- t_start:t_end
  phi_mat <- matrix(
    NA_real_,
    nrow = J,
    ncol = length(Tidx)
  )

  for (k_it in seq_along(Tidx)) {

    t = Tidx[k_it]
    ss  = ss_chain[, t]
    gs  = gs_chain[, t]
    cls <- as.integer(CL[, t]) + offset
    # number of examinees assigned to each latent class
    n_c <- tabulate(cls, nbins = C)
    # latent-class weights
    w_c <- switch(
      weight,
      uniform = rep(1 / C, C),
      pi = pi_chain[, t] / sum(pi_chain[, t]),
      nc_over_N = n_c / sum(n_c)
    )

    # compute item-specific dispersion for the current iteration
    for (j in 1:J) {

      # DINA conditional mean for each latent class
      eta_cj <- as.numeric(abs(AQ[, j] - qq[j]) < 1e-12)
      mu_cj <- ifelse(eta_cj == 1, 1 - ss[j], gs[j])
      mu_cj <- pmin(pmax(mu_cj, eps), 1 - eps)
      denom <- mu_cj * (1 - mu_cj)

      # observed proportion correct within each latent class
      ybar_cj <- numeric(C)
      for (c in 1:C) {
        idx <- which(cls == c)

        if (length(idx) > 0) {
          ybar_cj[c] <- mean(Y[idx, j])
        }
      }
      resid2 <- (ybar_cj - mu_cj)^2
      phi_mat[j, k_it] <- sum(
        w_c * (n_c * resid2 / denom)
      )
    }
  }

  # posterior summaries of item-specific dispersion
  phi_mean <- rowMeans(phi_mat, na.rm = TRUE)
  phi_median <- apply(phi_mat, 1, median, na.rm = TRUE)

  # posterior probability phi_j exceeds 1 + delta
  threshold <- 1 + delta
  p_phi_gt <- rowMeans(phi_mat > threshold, na.rm = TRUE)

  # select posterior summary and constrain dispersion to phi_j >= 1
  base_phi <- if (use_stat == "mean") phi_mean else phi_median
  phi_pos <- pmax(1, base_phi)

  # retain phi_j > 1 only if posterior evidence exceeds prob_cut
  cap_to_1 <- p_phi_gt < prob_cut
  phi_final <- phi_pos
  phi_final[cap_to_1] <- 1

  list(
    phi_by_iter = phi_mat,
    phi_mean = phi_mean,
    phi_median = phi_median,
    p_phi_gt = p_phi_gt,
    cap_to_1 = cap_to_1,
    phi_final = phi_final,
    prob_cut = prob_cut,
    delta = delta,
    use_stat = use_stat,
    weight = weight,
    t_start = t_start,
    t_end = t_end
  )
}

############################################################
# input files
############################################################
# the input is each folder with simulated datasets

path_rdatas  <- "./input/Simulation5/j12_n2000_pi50_A1clean30"
r_data_files <- list.files(path_rdatas, full.names = TRUE)

############################################################
# parallel configuration
############################################################
# I choose 15 or 20 workers

n_workers <- 15
cl <- makeCluster(n_workers)

clusterExport(
  cl,
  c("cpp_file", "r_data_files", "phi_BBQ_window"),
  envir = environment()
)

clusterEvalQ(cl, {
  library(Rcpp)
  library(RcppArmadillo)
  library(simcdm)
  library(dina)

  Rcpp::sourceCpp(cpp_file)

  NULL
})

############################################################
# Output directory
############################################################
# The output should coincide with the input. This can be 
# automataized, but my experiments were by hand.

out_dir <- "./output/Simulation5/j12_n2000_pi50_A1clean30"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clusterExport(cl, "out_dir", envir = environment())

############################################################
# Fit one simulated dataset
############################################################

run_one <- function(kk) {

  load(r_data_files[kk])  # Loads Y_sim, Q, and CLs

  # standard DINA model
  sim0 <- dina(
    Y = Y_sim,
    Q = Q,
    chain_length = 5000
  )

  # estimate item-specific dispersion parameters
  out_phi <- phi_BBQ_window(
    sim0,
    Y_sim,
    Q,
    t_start = 2000,
    t_end = 5000,
    delta = 0.5,
    weight = "uniform"
  )

  phi_bbq <- pmax(1, out_phi$phi_final)

  # qDINA model
  sim1 <- DINA_Gibbs_cpp_phi_item(
    Y = Y_sim,
    Q = Q,
    chain_length = 5000,
    phi = phi_bbq
  )

  # Save fitted models and calibration results
  tag <- tools::file_path_sans_ext(
    basename(r_data_files[kk])
  )

  save(
    CLs,
    sim0,
    sim1,
    out_phi,
    file = file.path(
      out_dir,
      paste0("fit_", tag, "_phiBBQ.RData")
    )
  )

  TRUE
}

############################################################
# Run simulation conditions in parallel
############################################################

tm <- system.time({
  ok <- parLapplyLB(
    cl,
    X = seq_along(r_data_files),
    fun = run_one
  )
})

stopCluster(cl)

tm
