set.seed(0)
library(rbiips)

#' # Compile model
modelfile <- 'inst/extdata/hmm.bug'
stopifnot(nchar(modelfile) > 0)
cat(readLines(modelfile), sep = '\n')

data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1), logtau = log(1))
model <- biips_model(modelfile, data, sample_data = TRUE)

#' # SMC algorithm
n_part <- 100
out_smc <- biips_smc_samples(model, c('x', 'c[2:10]'), n_part, type = 'fs',
                             rs_thres = 0.5, rs_type = 'stratified')

#' # PIMH algorithm
n_part <- 50
obj_pimh <- biips_pimh_init(model, c('x', 'c[2:10]'))  # Initialize
out_pimh_burn <- biips_pimh_update(obj_pimh, 100, n_part)  # Burn-in
out_pimh <- biips_pimh_samples(obj_pimh, 100, n_part)  # Samples

#' # PMMH algorithm
data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1))
model <- biips_model(modelfile, data)

obj_pmmh <- biips_pmmh_init(model, 'logtau', latent_names = c('x', 'c[2:10]'),
                            inits = list(logtau = -2))  # Initialize
out_pmmh_burn <- biips_pmmh_update(obj_pmmh, 100, n_part)  # Burn-in
out_pmmh <- biips_pmmh_samples(obj_pmmh, 100, n_part, thin = 1)  # Samples

#' # Save data
devtools::use_data(out_smc, out_pimh, out_pmmh, overwrite = TRUE)
