############################################################
# Simulation data generator
# qDINA simulation study
############################################################

# Libraries

library(simcdm)


############################################################
# data-generating mechanism
############################################################

sim_mixed_DINA_rebelA1 <- function(
  N, J, CLs, Q,
  g_good = 0.10,
  s_good = 0.10,
  g_A1clean = 0.20,
  s_A1clean = 0.20,
  bad_items,
  pi_010 = 0.30,
  pi_001 = 0.30,
  pi_011 = 0.30,
  p_rebel = 0.90
) {

  K <- ncol(Q)

  tbl_class <- simcdm::attribute_classes(K)
  A <- as.matrix(tbl_class[CLs + 1L, , drop = FALSE])

  Y <- matrix(0, N, J)

  # examinees with alpha_1 = 0 in the relevant latent strata
  R010 <- (A[, 1] == 0 & A[, 2] == 1)
  Z010 <- rbinom(N, 1, pi_010) * R010

  R001 <- (A[, 1] == 0 & A[, 3] == 1)
  Z001 <- rbinom(N, 1, pi_001) * R001

  R011 <- (A[, 1] == 0 & A[, 2] == 1 & A[, 3] == 1)
  Z011 <- rbinom(N, 1, pi_011) * R011

  for (j in seq_len(J)) {

    # Clean anchor item for attribute 2
    if (all(Q[j, ] == c(0, 1, 0))) {

      p <- ifelse(
        A[, 2] == 1,
        1 - s_good,
        g_good
      )

    # clean anchor item for attribute 3
    } else if (all(Q[j, ] == c(0, 0, 1))) {

      p <- ifelse(
        A[, 3] == 1,
        1 - s_good,
        g_good
      )

    # clean anchor item for attribute 1
    } else if (all(Q[j, ] == c(1, 0, 0))) {

      p <- ifelse(
        A[, 1] == 1,
        1 - s_A1clean,
        g_A1clean
      )

    # Misspecified items
    } else if (j %in% bad_items) {

      req <- which(Q[j, ] == 1)

      eta <- if (length(req) == 0) {
        rep(1, N)
      } else {
        apply(A[, req, drop = FALSE], 1, prod)
      }

      # baseline DINA response probability
      p <- ifelse(
        eta == 1,
        1 - s_good,
        g_good
      )

      # Misspecified 110 items:
      # selected examinees with pattern 010 behave as masters
      if (all(Q[j, ] == c(1, 1, 0))) {
        idx <- (eta == 0 & Z010 == 1)
        p[idx] <- p_rebel
      }

      # misspecified 101 items:
      # selected examinees with pattern 001 behave as masters
      if (all(Q[j, ] == c(1, 0, 1))) {
        idx <- (eta == 0 & Z001 == 1)
        p[idx] <- p_rebel
      }

      # Misspecified 111 item:
      # selected examinees with pattern 011 behave as masters
      if (all(Q[j, ] == c(1, 1, 1))) {
        idx <- (eta == 0 & Z011 == 1)
        p[idx] <- p_rebel
      }

    # Remaining items follow the standard DINA model
    } else {

      req <- which(Q[j, ] == 1)

      eta <- if (length(req) == 0) {
        rep(1, N)
      } else {
        apply(A[, req, drop = FALSE], 1, prod)
      }

      p <- ifelse(
        eta == 1,
        1 - s_good,
        g_good
      )
    }

    Y[, j] <- rbinom(N, 1, p)
  }

  Y
}


############################################################
# Generate all simulation conditions
############################################################

run_all_scenarios_sim5 <- function(
  base_dir = "./input/Simulation5",
  reps = 100,
  N_vec = c(500, 1000, 2000),
  pi_vec = c(0.30, 0.50),
  A1clean_vec = c(0.10, 0.30),
  J = 12,
  K = 3,
  g_good = 0.10,
  s_good = 0.10,
  p_rebel = 0.90,
  seed_base = 12345
) {

  # Q-matrix used in the simulation study
  Q <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1, 1, 0),
    c(1, 1, 0),
    c(1, 1, 0),
    c(1, 0, 1),
    c(1, 0, 1),
    c(1, 0, 1),
    c(0, 1, 1),
    c(0, 1, 1),
    c(1, 1, 1)
  )

  stopifnot(
    nrow(Q) == J,
    ncol(Q) == K
  )

  # Misspecified items
  bad_items <- c(4, 5, 7, 8, 12)

  # Simulation design
  scenarios <- expand.grid(
    N = N_vec,
    pi = pi_vec,
    A1clean = A1clean_vec,
    stringsAsFactors = FALSE
  )

  scenarios$tag <- sprintf(
    "j%d_n%d_pi%02d_A1clean%02d",
    J,
    scenarios$N,
    round(100 * scenarios$pi),
    round(100 * scenarios$A1clean)
  )

  scenarios$path <- file.path(
    base_dir,
    scenarios$tag
  )
  # Create one directory for each simulation condition
  for (p in scenarios$path) {
    dir.create(
      p,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  # Generate replications for each condition
  for (s in seq_len(nrow(scenarios))) {

    N <- scenarios$N[s]
    pi <- scenarios$pi[s]
    A1clean <- scenarios$A1clean[s]
    pathw <- scenarios$path[s]

    # Uniform true latent-class distribution
    PIs_true <- rep(
      1 / (2^K),
      2^K
    )

    for (jj in seq_len(reps)) {

      # Condition and Replication specific random seed
      set.seed(
        seed_base + 100000 * s + jj
      )

      CLs <- sample.int(
        2^K,
        size = N,
        replace = TRUE,
        prob = PIs_true
      ) - 1L
      CLs <- as.integer(CLs)

      Y_sim <- sim_mixed_DINA_rebelA1(
        N = N,
        J = J,
        CLs = CLs,
        Q = Q,
        g_good = g_good,
        s_good = s_good,
        g_A1clean = A1clean,
        s_A1clean = A1clean,
        bad_items = bad_items,
        pi_010 = pi,
        pi_001 = pi,
        pi_011 = pi,
        p_rebel = p_rebel
      )
      fout <- file.path(
        pathw,
        sprintf(
          "s5_%s_rep_%03d.RData",
          scenarios$tag[s],
          jj - 1L
        )
      )

      save(
        Y_sim,
        Q,
        CLs,
        bad_items,
        g_good,
        s_good,
        A1clean,
        pi,
        p_rebel,
        file = fout
      )
    }
  }

  invisible(scenarios)
}


############################################################
# Run simulation data generation
############################################################

scenarios_created <- run_all_scenarios_sim5(
  base_dir = "./input/Simulation5",
  reps = 100,
  N_vec = c(500, 1000, 2000),
  pi_vec = c(0.30, 0.50),
  A1clean_vec = c(0.10, 0.30),
  J = 12,
  K = 3,
  g_good = 0.10,
  s_good = 0.10,
  p_rebel = 0.90
)

scenarios_created[, c("tag", "path")]
