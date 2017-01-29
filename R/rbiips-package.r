release_questions <- function() {
  c(
    "Have you updated the copyright year in DESCRIPTION?",
    "Have you updated the docs? devtools::document()",
    "Have you updated the website? pkgdown::build_site()",
    "Have you pulled biips submodule? cd src/biips; git pull"
  )
}

#' Bayesian inference with interacting particle systems
#'
#' \pkg{rbiips} is an interface with the \href{https://biips.github.io}{Biips} C++ libraries for analysing
#' Bayesian graphical models using advanced particle methods.
#'
#' Biips is a general software for Bayesian inference with interacting particle
#' systems, a.k.a. sequential Monte Carlo (SMC) methods. It aims at popularizing
#' the use of these methods to non-statistician researchers and students, thanks
#' to its automated \dQuote{black box} inference engine. It borrows from the
#' \href{http://www.mrc-bsu.cam.ac.uk/software/bugs/}{BUGS}/\href{http://mcmc-jags.sourceforge.net/}{JAGS}
#' software, widely used in Bayesian statistics, the statistical
#' modeling with graphical models and the language associated with their
#' descriptions.
#'
#' See the \href{https://biips.github.io}{Biips website} for more
#' information.
#'
#' The typical workflow is the following:
#' \itemize{
#'   \item Define the model in BUGS language (see the \href{http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}{JAGS User Manual}
#'   for help) and the data.
#'   \item Add custom functions or distributions with
#'     \code{\link{biips_add_function}} and \code{\link{biips_add_distribution}}.
#'   \item Compile the model with \code{\link{biips_model}}
#'   \item Run inference algorithms:
#'     \itemize{
#'       \item Analyse sensitivity to parameters with \code{\link{biips_smc_sensitivity}}.
#'       \item Run SMC filtering and smoothing algorithms with \code{\link{biips_smc_samples}}.
#'       \item Run particle MCMC algorithms with \code{\link{biips_pimh_samples}} or
#'         \code{\link{biips_pmmh_samples}}.
#'     }
#'   \item Diagnose and analyze the output obtained as \code{\link{smcarray}} and
#'     \code{\link{mcmcarray}} objects with \code{\link{biips_diagnosis}},
#'     \code{\link{biips_summary}}, \code{\link{biips_density}},
#'     \code{\link{biips_hist}} and \code{\link{biips_table}}
#'  }
#'
#' @name rbiips-package
#' @aliases rbiips
#' @docType package
#' @author 
#'   \itemize{
#'     \item \href{http://adrien.tspace.fr}{Adrien
#'     Todeschini}
#'     \item \href{http://www.stats.ox.ac.uk/~caron/}{François Caron}
#'     \item Marc Fuentes
#'   }
#'   \pkg{rbiips} development was initiated by the research team
#'   \href{http://alea.bordeaux.inria.fr}{ALEA} at
#'   \href{http://www.inria.fr/en/centre/bordeaux}{Inria Bordeaux Sud-Ouest}.
#' @seealso \code{\link{biips_add_function}}, \code{\link{biips_add_distribution}},
#'   \code{\link{biips_model}}, \code{\link{biips_smc_sensitivity}}, \code{\link{biips_smc_samples}},
#'   \code{\link{biips_pimh_samples}}, \code{\link{biips_pmmh_samples}}, \code{\link{smcarray}},
#'   \code{\link{mcmcarray}}, \code{\link{biips_diagnosis}}, \code{\link{biips_summary}},
#'   \code{\link{biips_density}}, \code{\link{biips_hist}}, \code{\link{biips_table}},
#'   \href{https://biips.github.io}{Biips website},
#'   \href{http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}{JAGS User Manual}
#' @references Adrien Todeschini, François Caron, Marc Fuentes, Pierrick Legrand, Pierre Del Moral (2014).
#'   Biips: Software for Bayesian Inference with Interacting Particle Systems. 
#'   arXiv preprint arXiv:1412.3779.
#'   URL \url{http://arxiv.org/abs/1412.3779}
#' @keywords package
#' @useDynLib rbiips
#' @importFrom Rcpp evalCpp
#' @examples
#' #' # Add custom functions and distributions to BUGS language
#' #' Add custom function `f`
#' f_dim <- function(x_dim, t_dim) {
#'   # Check dimensions of the input and return dimension of the output of function f
#'   stopifnot(prod(x_dim) == 1, prod(t_dim) == 1)
#'   x_dim
#' }
#' f_eval <- function(x, t) {
#'   # Evaluate function f
#'   0.5 * x + 25 * x/(1 + x^2) + 8 * cos(1.2 * t)
#' }
#' biips_add_function('f', 2, f_dim, f_eval)
#'
#' #' Add custom sampling distribution `dMN`
#' dMN_dim <- function(mu_dim, Sig_dim) {
#'   # Check dimensions of the input and return dimension of the output of
#'   # distribution dMN
#'   stopifnot(prod(mu_dim) == mu_dim[1], length(Sig_dim) == 2, mu_dim[1] == Sig_dim)
#'   mu_dim
#' }
#' dMN_sample <- function(mu, Sig) {
#'   # Draw a sample of distribution dMN
#'   mu + t(chol(Sig)) %*% rnorm(length(mu))
#' }
#' biips_add_distribution('dMN', 2, dMN_dim, dMN_sample)
#'
#' #' # Compile model
#' modelfile <- system.file('extdata', 'hmm_f.bug', package = 'rbiips')
#' stopifnot(nchar(modelfile) > 0)
#' cat(readLines(modelfile), sep = '\n')
#'
#' data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1), logtau = log(1))
#' model <- biips_model(modelfile, data, sample_data = TRUE)
#'
#' #' # SMC algorithm
#' n_part <- 100
#' out_smc <- biips_smc_samples(model, c('x', 'c[2:10]'), n_part, type = 'fs',
#'                              rs_thres = 0.5, rs_type = 'stratified')
#'
#' biips_diagnosis(out_smc)
#' biips_summary(out_smc)
#' par(mfrow = c(2, 2))
#' plot(biips_density(out_smc$x, bw = 'nrd0', adjust = 1, n = 100))
#' plot(biips_table(out_smc[['c[2:10]']]))
"_PACKAGE"
