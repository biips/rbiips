rbiips
========
[![Travis-CI Linux](https://travis-matrix-badges.herokuapp.com/repos/biips/rbiips/branches/master/1)](https://travis-ci.org/biips/rbiips)
[![Travis-CI OS X](https://travis-matrix-badges.herokuapp.com/repos/biips/rbiips/branches/master/2)](https://travis-ci.org/biips/rbiips)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/biips/rbiips?branch=master&svg=true)](https://ci.appveyor.com/project/biips/rbiips)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rbiips)](https://cran.r-project.org/package=rbiips)
[![GPLv3 License](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

Bayesian inference with interacting particle systems

**rbiips** is an interface with the [Biips][Biips website] C++ libraries for analysing
Bayesian graphical models using advanced particle methods.

Biips is a general software for Bayesian inference with interacting particle
systems, a.k.a. sequential Monte Carlo (SMC) methods. It aims at popularizing
the use of these methods to non-statistician researchers and students, thanks
to its automated "black box" inference engine. It borrows from the
[BUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/)/[JAGS][JAGS website]
software, widely used in Bayesian statistics, the statistical
modeling with graphical models and the language associated with their
descriptions.

The typical workflow is the following:

- Define the model in BUGS language (see the [JAGS User Manual][JAGS User Manual]
  for help) and the data.
- Add custom functions or distributions with 
  [`biips_add_function`](https://biips.github.io/rbiips/reference/extend-bugs.html) and 
  [`biips_add_distribution`](https://biips.github.io/rbiips/reference/extend-bugs.html).
- Compile the model with [`biips_model`](https://biips.github.io/rbiips/reference/biips_model.html)
- Run inference algorithms:
    - Analyse sensitivity to parameters with 
      [`biips_smc_sensitivity`](https://biips.github.io/rbiips/reference/biips_smc_sensitivity.html).
    - Run SMC filtering and smoothing algorithms with 
      [`biips_smc_samples`](https://biips.github.io/rbiips/reference/biips_smc_samples.html).
    - Run particle MCMC algorithms with 
      [`biips_pimh_samples`](https://biips.github.io/rbiips/reference/pimh-object.html) or
      [`biips_pmmh_samples`](https://biips.github.io/rbiips/reference/pmmh-object.html).
- Diagnose and analyze the output obtained as 
  [`smcarray`](https://biips.github.io/rbiips/reference/smcarray-object.html) and 
  [`mcmcarray`](https://biips.github.io/rbiips/reference/mcmcarray-object.html) objects with 
  [`biips_summary`](https://biips.github.io/rbiips/reference/smcarray-object.html), 
  [`biips_density`](https://biips.github.io/rbiips/reference/smcarray-object.html),
  [`biips_hist`](https://biips.github.io/rbiips/reference/smcarray-object.html) and 
  [`biips_table`](https://biips.github.io/rbiips/reference/smcarray-object.html)

Installation
------------
Install the latest version from [GitHub](https://github.com/biips/rbiips):
```R
install_dir <- tempdir()
system(paste("git clone --recursive", shQuote("https://github.com/biips/rbiips.git"), shQuote(install_dir)))
devtools::install(install_dir)
```

Authors
----------

- [Adrien Todeschini](http://adrien.tspace.fr)
- [François Caron](http://www.stats.ox.ac.uk/~caron/)
- Marc Fuentes

**rbiips** development was initiated by the research team
[ALEA](http://alea.bordeaux.inria.fr) at
[Inria Bordeaux Sud-Ouest](http://www.inria.fr/en/centre/bordeaux).

References
-------------
Adrien Todeschini, François Caron, Marc Fuentes, Pierrick Legrand, Pierre Del Moral (2014). 
Biips: Software for Bayesian Inference with Interacting Particle Systems. 
arXiv preprint arXiv:1412.3779. 
URL <http://arxiv.org/abs/1412.3779>.

[Biips website]: https://biips.github.io
[JAGS website]: http://mcmc-jags.sourceforge.net/
[JAGS User Manual]: http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download

Example 
----------
```R
library(rbiips)

#' # Add custom functions and distributions to BUGS language
#' Add custom function `f`
f_dim <- function(x_dim, t_dim) {
  # Check dimensions of the input and return dimension of the output of function f
  stopifnot(prod(x_dim) == 1, prod(t_dim) == 1)
  x_dim
}
f_eval <- function(x, t) {
  # Evaluate function f
  0.5 * x + 25 * x/(1 + x^2) + 8 * cos(1.2 * t)
}
biips_add_function('f', 2, f_dim, f_eval)

#' Add custom sampling distribution `dMN`
dMN_dim <- function(mu_dim, Sig_dim) {
  # Check dimensions of the input and return dimension of the output of
  # distribution dMN
  stopifnot(prod(mu_dim) == mu_dim[1], length(Sig_dim) == 2, mu_dim[1] == Sig_dim)
  mu_dim
}
dMN_sample <- function(mu, Sig) {
  # Draw a sample of distribution dMN
  mu + t(chol(Sig)) %*% rnorm(length(mu))
}
biips_add_distribution('dMN', 2, dMN_dim, dMN_sample)

#' # Compile model
modelfile <- system.file('extdata', 'hmm_f.bug', package = 'rbiips')
stopifnot(nchar(modelfile) > 0)
cat(readLines(modelfile), sep = '\n')

data_ <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1), logtau = log(1))
model <- biips_model(modelfile, data_, sample_data = TRUE)

#' # SMC algorithm
n_part <- 100
out_smc <- biips_smc_samples(model, c('x', 'c[2:10]'), n_part, type = 'fs',
                             rs_thres = 0.5, rs_type = 'stratified')

biips_diagnosis(out_smc)
biips_summary(out_smc)
par(mfrow = c(2, 2))
plot(biips_density(out_smc$x, bw = 'nrd0', adjust = 1, n = 100))
plot(biips_table(out_smc[['c[2:10]']]))

#' # PIMH algorithm
n_part <- 50
obj_pimh <- biips_pimh_init(model, c('x', 'c[2:10]'))  # Initialize
out_pimh_burn <- biips_pimh_update(obj_pimh, 100, n_part)  # Burn-in
out_pimh <- biips_pimh_samples(obj_pimh, 100, n_part)  # Samples

biips_summary(out_pimh)
par(mfrow = c(2, 2))
plot(biips_density(out_pimh$x))
biips_hist(out_pimh$x)
plot(biips_table(out_pimh[['c[2:10]']]))

#' # SMC sensitivity analysis
n_part <- 50
logtau_val <- -10:10
out_sens <- biips_smc_sensitivity(model, list(logtau = logtau_val), n_part)

#' # PMMH algorithm
data_ <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1))
model <- biips_model(modelfile, data_)

n_part <- 50
obj_pmmh <- biips_pmmh_init(model, 'logtau', latent_names = c('x', 'c[2:10]'),
                            inits = list(logtau = -2))  # Initialize
out_pmmh_burn <- biips_pmmh_update(obj_pmmh, 100, n_part)  # Burn-in
out_pmmh <- biips_pmmh_samples(obj_pmmh, 100, n_part, thin = 1)  # Samples

biips_summary(out_pmmh)
par(mfrow = c(2, 2))
plot(biips_density(out_pmmh$logtau))
biips_hist(out_pmmh$logtau)
plot(biips_density(out_pmmh$x))
biips_hist(out_pmmh$x)
plot(biips_table(out_pmmh[['c[2:10]']]))
```
