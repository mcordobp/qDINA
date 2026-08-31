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
source("./R/phi_calibration.R")
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