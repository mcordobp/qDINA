// -----------------------------------------------------------------------------
// qDINA Gibbs sampler
//
// This implementation is adapted from the Gibbs sampling framework for the
// DINA model described by Culpepper (2015) and its implementation in the
// R package "dina" by Culpepper and Balamuta. The sampler is extended here to
// incorporate item-specific dispersion parameters, phi_j, through power
// weights w_j = 1 / phi_j.
//
// References:
// Culpepper, S. A. (2015). Bayesian estimation of the DINA model with Gibbs
// sampling. Journal of Educational and Behavioral Statistics, 40(5), 454–476.
// https://doi.org/10.3102/1076998615595403
//
// Culpepper, S. A., & Balamuta, J. J. (2019). dina: Bayesian estimation of
// DINA model. R package. https://cran.r-project.org/package=dina
// -----------------------------------------------------------------------------

#include <cmath>
#include <RcppArmadillo.h>
#include <simcdm.h>
#include <rgen.h>

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List update_alpha_phi(const arma::mat &Amat, const arma::mat &Q,
                            const arma::vec &ss, const arma::vec &gs,
                            const arma::mat &Y, const arma::vec &PIs,
                            arma::mat &ALPHAS, const arma::vec &delta0,
                            const arma::vec &phi) {

  const unsigned int N = Y.n_rows;
  const unsigned int J = Y.n_cols;
  const unsigned int C = Amat.n_rows;
  const unsigned int K = Q.n_cols;

  if (phi.n_elem != J) Rcpp::stop("phi must have length J (number of items).");
  if (arma::any(phi <= 0.0)) Rcpp::stop("All phi_j must be > 0.");

  const double eps = 1e-12;

  arma::vec logPY(C);
  arma::vec PS(C);
  arma::vec Ncs = arma::zeros<arma::vec>(C);
  arma::vec CLASSES(N);

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int c = 0; c < C; c++) {

      double llw = 0.0; // Log-likelihood weighted by 1/phi_j

      for (unsigned int j = 0; j < J; j++) {

        // eta_ij = 1 if the latent class masters all attributes required by item j
        double etaij = 1.0;
        if (arma::dot(Amat.row(c), Q.row(j)) < arma::dot(Q.row(j), Q.row(j))) {
          etaij = 0.0;
        }

        const double y = Y(i, j);
        double p;

        if (etaij == 1.0 && y == 1.0)       p = 1.0 - ss(j);
        else if (etaij == 0.0 && y == 1.0)  p = gs(j);
        else if (etaij == 1.0 && y == 0.0)  p = ss(j);
        else                                p = 1.0 - gs(j);

        p = std::min(std::max(p, eps), 1.0 - eps);

        llw += (1.0 / phi(j)) * std::log(p);
      }

      logPY(c) = std::log(std::max(PIs(c), eps)) + llw;
    }

    // Numerically stable softmax
    double m = logPY.max();
    arma::vec tmp = arma::exp(logPY - m);
    PS = tmp / arma::accu(tmp);

    double ci = rgen::rmultinomial(PS);
    ALPHAS.row(i) = Amat.row((unsigned)ci);
    Ncs((unsigned)ci) += 1.0;
    CLASSES(i) = ci;
  }

  arma::vec PIs_new = rgen::rdirichlet(Ncs + delta0);

  return Rcpp::List::create(
    Rcpp::Named("PS") = PS,
    Rcpp::Named("PIs_new") = PIs_new,
    Rcpp::Named("CLASSES") = CLASSES
  );
}


// [[Rcpp::export]]
Rcpp::List update_sg_phi(const arma::mat &Y, const arma::mat &Q,
                         const arma::mat &ALPHAS, const arma::vec &ss_old,
                         double as0, double bs0, double ag0, double bg0,
                         const arma::vec &phi) {

  const unsigned int N = Y.n_rows;
  const unsigned int J = Y.n_cols;

  if (phi.n_elem != J) Rcpp::stop("phi must have length J (number of items).");
  if (arma::any(phi <= 0.0)) Rcpp::stop("All phi_j must be > 0.");

  const double eps = 1e-12;

  arma::vec ss_new(J);
  arma::vec gs_new(J);

  arma::mat AQ = ALPHAS * Q.t();

  for (unsigned int j = 0; j < J; j++) {
    const double wj = 1.0 / phi(j);

    double us = R::runif(0, 1);
    double ug = R::runif(0, 1);

    arma::vec ETA = arma::zeros<arma::vec>(N);
    double qq = arma::as_scalar(Q.row(j) * (Q.row(j)).t());
    ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);

    double y_dot_eta = arma::as_scalar((Y.col(j)).t() * ETA);

    double T = arma::accu(ETA);                  // Masters all required attributes
    double S = T - y_dot_eta;                    // Slips: eta = 1, y = 0
    double G = arma::accu(Y.col(j)) - y_dot_eta; // Guesses: eta = 0, y = 1

    double D = (double)N - T - G; // eta = 0, y = 0
    double B = y_dot_eta;         // eta = 1, y = 1

    // Tempered Beta parameters
    double a_g = ag0 + wj * G;
    double b_g = bg0 + wj * D;

    double a_s = as0 + wj * S;
    double b_s = bs0 + wj * B;

    // g | s_old truncated to (0, 1 - s_old)
    double upper_g = 1.0 - ss_old(j);
    upper_g = std::min(std::max(upper_g, eps), 1.0 - eps);
    double pg = R::pbeta(upper_g, a_g, b_g, 1, 0);
    pg = std::min(std::max(pg, eps), 1.0 - eps);
    gs_new(j) = R::qbeta(ug * pg, a_g, b_g, 1, 0);

    // s | g_new truncated to (0, 1 - g_new)
    double upper_s = 1.0 - gs_new(j);
    upper_s = std::min(std::max(upper_s, eps), 1.0 - eps);
    double ps = R::pbeta(upper_s, a_s, b_s, 1, 0);
    ps = std::min(std::max(ps, eps), 1.0 - eps);
    ss_new(j) = R::qbeta(us * ps, a_s, b_s, 1, 0);
  }

  return Rcpp::List::create(
    Rcpp::Named("ss_new") = ss_new,
    Rcpp::Named("gs_new") = gs_new
  );
}


// [[Rcpp::export]]
Rcpp::List DINA_Gibbs_cpp_phi_item(const arma::mat &Y, const arma::mat &Q,
                                   const arma::vec &phi,
                                   unsigned int chain_length = 10000) {

  const unsigned int N = Y.n_rows;
  const unsigned int J = Y.n_cols;
  const unsigned int K = Q.n_cols;
  const unsigned int C = (unsigned int)std::pow(2.0, (double)K);

  if (phi.n_elem != J) Rcpp::stop("phi must have length J (number of items).");
  if (arma::any(phi <= 0.0)) Rcpp::stop("All phi_j must be > 0.");

  arma::mat Amat = simcdm::attribute_classes(K);

  arma::vec delta0 = arma::ones<arma::vec>(C);
  double as0 = 1.0, bs0 = 1.0, ag0 = 1.0, bg0 = 1.0;

  arma::mat SigS(J, chain_length);
  arma::mat GamS(J, chain_length);
  arma::mat PIs(C, chain_length);
  arma::mat CLASSES(N, chain_length);

  arma::mat alphas = arma::randu<arma::mat>(N, K);
  alphas.elem(find(alphas > 0.5)).ones();
  alphas.elem(find(alphas <= 0.5)).zeros();

  arma::vec ss = arma::randu<arma::vec>(J);
  arma::vec gs = arma::randu<arma::vec>(J);
  arma::vec pis = arma::randu<arma::vec>(C);

  for (unsigned int t = 0; t < chain_length; t++) {

    // Update latent classes and class proportions
    Rcpp::List step1a =
      update_alpha_phi(Amat, Q, ss, gs, Y, pis, alphas, delta0, phi);

    pis = Rcpp::as<arma::vec>(step1a["PIs_new"]);
    CLASSES.col(t) = Rcpp::as<arma::vec>(step1a["CLASSES"]);
    PIs.col(t) = pis;

    // Update slipping and guessing parameters
    Rcpp::List step1b =
      update_sg_phi(Y, Q, alphas, ss, as0, bs0, ag0, bg0, phi);

    ss = Rcpp::as<arma::vec>(step1b["ss_new"]);
    gs = Rcpp::as<arma::vec>(step1b["gs_new"]);

    SigS.col(t) = ss;
    GamS.col(t) = gs;
  }

  return Rcpp::List::create(
    Rcpp::Named("CLASSES") = CLASSES,
    Rcpp::Named("PIs") = PIs,
    Rcpp::Named("SigS") = SigS,
    Rcpp::Named("GamS") = GamS,
    Rcpp::Named("phi") = phi
  );
}
