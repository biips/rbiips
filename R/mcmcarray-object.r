#' @title Objects for representing MCMC output.
#'
#' @description A \code{mcmcarray} object is returned by the
#'   \code{\link{biips_pimh_samples}} or \code{\link{biips_pmmh_samples}}
#'   functions to represent MCMC output of a given variable.
#'
#'   A \code{mcmcarray.list} object is a named list of \code{mcmcarray} objects
#'   for different monitored variables.
#'
#'   The methods apply identically to \code{mcmcarray} or \code{mcmcarray.list}
#'   objects and return a named list with the same named members as the input
#'   object.
#'
#' @name mcmcarray-object
#' @aliases mcmcarray.list-object mcmcarray mcmcarray.list
#' @param object,x a \code{mcmcarray} or \code{mcmcarray.list} object.
#' @param ... additional arguments to be passed to the default methods. See
#'   \code{\link[stats]{density}}, \code{\link[graphics]{hist}},
#'   \code{\link[base]{table}}
#'
#' @return The methods apply identically to \code{mcmcarray} or
#'   \code{mcmcarray.list} objects and return a named list with the same named
#'   members as the input object.
#'
#' @examples
#' modelfile <- system.file('extdata', 'hmm.bug', package = 'rbiips')
#' stopifnot(nchar(modelfile) > 0)
#' cat(readLines(modelfile), sep = '\n')
#'
#' #' # PIMH algorithm
#' data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1), logtau = log(1))
#' model <- biips_model(modelfile, data, sample_data = TRUE)
#'
#' n_part <- 50
#' obj_pimh <- biips_pimh_init(model, c('x', 'c[2:10]'))  # Initialize
#' out_pimh_burn <- biips_pimh_update(obj_pimh, 100, n_part)  # Burn-in
#' out_pimh <- biips_pimh_samples(obj_pimh, 100, n_part)  # Samples
#'
#' #' Manipulate `mcmcarray.list` object
#' is.mcmcarray.list(out_pimh)
#' names(out_pimh)
#' out_pimh
#' biips_summary(out_pimh)
#'
#' #' Manipulate `mcmcarray` object
#' is.mcmcarray(out_pimh$x)
#' out_pimh$x
#' summ_pimh_x <- biips_summary(out_pimh$x, order = 2, probs = c(0.025, 0.975))
#' summ_pimh_x
#' dens_pimh_x <- biips_density(out_pimh$x)
#' par(mfrow = c(2, 2))
#' plot(dens_pimh_x)
#' par(mfrow = c(2, 2))
#' biips_hist(out_pimh$x)
#'
#' is.mcmcarray(out_pimh[['c[2:10]']])
#' out_pimh[['c[2:10]']]
#' summ_pimh_c <- biips_summary(out_pimh[['c[2:10]']])
#' summ_pimh_c
#' table_pimh_c <- biips_table(out_pimh[['c[2:10]']])
#' par(mfrow = c(2, 2))
#' plot(table_pimh_c)
#'
#' #' # PMMH algorithm
#' data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1))
#' model <- biips_model(modelfile, data)
#'
#' n_part <- 50
#' obj_pmmh <- biips_pmmh_init(model, 'logtau', latent_names = c('x', 'c[2:10]'),
#'                             inits = list(logtau = -2))  # Initialize
#' out_pmmh_burn <- biips_pmmh_update(obj_pmmh, 100, n_part)  # Burn-in
#' out_pmmh <- biips_pmmh_samples(obj_pmmh, 100, n_part, thin = 1)  # Samples
#'
#' #' Manipulate `mcmcarray.list` object
#' is.mcmcarray.list(out_pmmh)
#' names(out_pmmh)
#' out_pmmh
#' biips_summary(out_pmmh)
#'
#' #' Manipulate `mcmcarray` object
#' is.mcmcarray(out_pmmh$logtau)
#' out_pmmh$logtau
#' summ_pmmh_lt <- biips_summary(out_pmmh$logtau, order = 2, probs = c(0.025, 0.975))
#' dens_pmmh_lt <- biips_density(out_pmmh$logtau)
#' par(mfrow = c(2, 1))
#' plot(dens_pmmh_lt)
#' biips_hist(out_pmmh$logtau)
#'
#' is.mcmcarray(out_pmmh$x)
#' out_pmmh$x
#' summ_pmmh_x <- biips_summary(out_pmmh$x, order = 2, probs = c(0.025, 0.975))
#' dens_pmmh_x <- biips_density(out_pmmh$x)
#' par(mfrow = c(2, 2))
#' plot(dens_pmmh_x)
#' par(mfrow = c(2, 2))
#' biips_hist(out_pmmh$x)
#'
#' is.mcmcarray(out_pmmh[['c[2:10]']])
#' out_pmmh[['c[2:10]']]
#' summ_pmmh_c <- biips_summary(out_pmmh[['c[2:10]']])
#' table_pmmh_c <- biips_table(out_pmmh[['c[2:10]']])
#' par(mfrow = c(2, 2))
#' plot(table_pmmh_c)
NULL

#' @rdname mcmcarray-object
#' @param data      numerical vector
#' @param dim       vector of integers. dimension of the array
#' @param dimnames  character vector
#' @param iteration integer. index of the dimension corresponding to iterations
#'   of the MCMC.
#' @param chain     integer. index of the dimension corresponding to chain of
#'   the MCMC.
#' @param name   string. variable name
#' @param lower   vector of integers. variable lower bound
#' @param upper   vector of integers. variable upper bound
#' @return The \code{mcmcarray} function returns an object of class \code{mcmcarray}.
mcmcarray <- function(data = NA, dim = length(data), dimnames = NULL, iteration = length(dim),
  chain = NA, name = "mcmcarray", lower = NULL, upper = NULL) {
  stopifnot(length(iteration) == 1, length(chain) == 1)
  if (is.na(iteration))
    stopifnot(is.na(chain)) else {
    stopifnot(is.numeric(iteration), iteration >= 1, iteration <= length(dim))
    if (!any(is.na(chain)))
      stopifnot(is.numeric(chain), chain >= 1, chain <= length(dim), chain !=
        iteration)
  }

  x <- array(data, dim, dimnames)

  nd <- names(dim(x))
  nd[iteration] <- "iteration"
  nd[chain] <- "chain"
  names(dim(x)) <- nd

  drop_dims <- names(dim(x)) %in% c("iteration", "chain")
  dim_array <- dim(x)[!drop_dims]

  if (is.null(lower)) {
    lower <- rep(1, length(dim_array))
  }
  if (is.null(upper)) {
    upper <- dim_array
  }

  attr(x, "name") <- name
  attr(x, "lower") <- lower
  attr(x, "upper") <- upper

  class(x) <- "mcmcarray"
  return(x)
}


#' @export
#' @rdname mcmcarray-object
#' @return The function \code{is.mcmcarray} returns \code{TRUE} if the object is
#'   of class \code{mcmcarray}.
is.mcmcarray <- function(object) {
  return(class(object) == "mcmcarray")
}

#' @export
#' @rdname mcmcarray-object
#' @return The function \code{is.mcmcarray.list} returns \code{TRUE} if the
#'   object is of class \code{mcmcarray.list}.
is.mcmcarray.list <- function(object) {
  return(class(object) == "mcmcarray.list")
}

#' @export
print.mcmcarray <- function(x, ...) {
  print(summary(x), ...)
}

#' @export
print.mcmcarray.list <- function(x, ...) {
  print(summary(x), ...)
}

#' @export
#' @rdname mcmcarray-object
#' @inheritParams biips_summary.smcarray
#' @return The method \code{biips_summary} returns univariate marginal
#'   statistics. The output innermost members are objects of class
#'   \code{summary.mcmcarray}, \emph{i.e.} lists with members:
#'   \item{mean}{mean, if \code{order>=1}.}
#'   \item{var}{variance, if \code{order>=2}.}
#'   \item{skew}{skewness, if \code{order>=3}.}
#'   \item{kurt}{kurtosis, if \code{order>=4}.}
#'   \item{probs}{vector of quantile probabilities.}
#'   \item{quant}{list of quantile values, if \code{probs} is not empty.}
#'   \item{mode}{most frequent values for discrete components.}
#' @importFrom stats var quantile
biips_summary.mcmcarray <- function(object, probs = c(), order = ifelse(mode, 0,
  1), mode = all(object == as.integer(object)), ...) {
  stopifnot(is.mcmcarray(object))

  ### TODO check arguments
  if (length(probs) > 0)
    stopifnot(is.numeric(probs), probs > 0, probs < 1)

  drop_dims <- names(dim(object)) %in% c("iteration", "chain")
  dim_array <- dim(object)[!drop_dims]
  marg <- which(!drop_dims)
  summ <- list()

  ### moment based statistics
  if (order >= 1)
    summ$mean <- apply(object, marg, FUN = mean)

  if (order >= 2)
    summ$var <- apply(object, marg, FUN = var)

  if (order >= 3) {
    mom2 <- apply(object^2, marg, FUN = mean)
    mom3 <- apply(object^3, marg, FUN = mean)
    summ$skew <- (mom3 - 3 * mom2 * summ$mean + 2 * summ$mean^3)/summ$var^(3/2)
  }

  if (order >= 4) {
    mom4 <- apply(object^4, marg, FUN = mean)
    summ$kurt <- (mom4 - 4 * mom3 * summ$mean + 6 * mom2 * summ$mean^2 - 3 *
      summ$mean^4)/summ$var^2 - 3
  }

  ### quantiles
  if (length(probs) > 0) {
    summ$probs <- probs
    quant <- apply(object, marg, FUN = quantile, probs = probs, names = FALSE)
    ndim <- length(dim(quant))
    quant <- aperm(quant, c(2:ndim, 1))

    len <- prod(dim_array)
    summ$quant <- list()
    for (i in 1:length(probs)) {
      summ$quant[[as.character(probs[i])]] <- array(quant[((i - 1) * len +
        1):(i * len)], dim = dim_array)
    }
  }

  ### mode
  if (mode) {
    summ$mode <- apply(object, marg, FUN = function(x) as.numeric(names(which.max(table(x)))))
  }

  summ$drop.dims <- dim(object)[drop_dims]

  class(summ) <- "summary.mcmcarray"

  return(summ)
}

#' @export
#' @rdname mcmcarray-object
biips_summary.mcmcarray.list <- function(object, ...) {
  stopifnot(is.mcmcarray.list(object))
  ans <- list()
  for (n in names(object)) {
    if (!is.mcmcarray(object[[n]]))
      next
    ans[[n]] <- biips_summary(object[[n]], ...)
  }

  class(ans) <- "summary.mcmcarray.list"

  return(ans)
}

#' @export
print.summary.mcmcarray <- function(x, ...) {
  cat("mcmcarray:\n")
  print(x[!(names(x) %in% c("drop.dims"))], ...)
  if (length(x$drop.dims) > 0) {
    cat("Marginalizing over:", paste(paste(names(x$drop.dims), "(", x$drop.dims,
      ")", sep = ""), collapse = ","), "\n")
  }
}


#' @export
print.summary.mcmcarray.list <- function(x, ...) {
  for (n in names(x)) {
    if (n %in% c("log_marg_like_pen", "log_marg_like", "info"))
      next
    cat(n, " ", sep = "")
    print(x[[n]], ...)
    cat("\n")
  }
  invisible()
}

#' @export
#' @rdname mcmcarray-object
#' @return The method \code{biips_table} returns univariate marginal frequency
#'   tables or probability mass estimates of discrete variables. The output
#'   innermost members are objects of class \code{table.mcmcarray}.
biips_table.mcmcarray <- function(x, ...) {
  stopifnot(is.mcmcarray(x))
  out <- list()

  dimen <- dim(x)
  drop_dim <- names(dimen) %in% c("iteration", "chain")
  n_samples <- prod(dimen[drop_dim])
  len <- prod(dimen[!drop_dim])

  for (d in 1:len) {
    ind_vec <- seq(d, len * (n_samples - 1) + d, len)
    values <- x[ind_vec]

    l <- attr(x, "lower")
    u <- attr(x, "upper")
    ind <- l + get_index(d, rep(0, length(l)), u - l)
    dnn <- deparse_varname(attr(x, "name"), ind, ind)

    out[[d]] <- table(values, dnn = dnn, ...)/n_samples
  }

  dim(out) <- dimen[!drop_dim]
  class(out) <- "table.mcmcarray"
  return(out)
}


#' @export
#' @rdname mcmcarray-object
#' @inheritParams biips_density.smcarray
#' @return The method \code{biips_density} returns univariate marginal kernel
#'   density estimates. The output innermost members are objects of class
#'   \code{density.mcmcarray}.
biips_density.mcmcarray <- function(x, bw = "nrd0", ...) {
  stopifnot(is.mcmcarray(x))
  out <- list()

  dimen <- dim(x)
  drop_dim <- names(dimen) %in% c("iteration", "chain")
  n_samples <- prod(dimen[drop_dim])
  len <- prod(dimen[!drop_dim])

  for (d in 1:len) {
    ind_vec <- seq(d, len * (n_samples - 1) + d, len)
    values <- x[ind_vec]

    out[[d]] <- density(values, bw = rec(bw, d), ...)  # recycle bw
    ind <- attr(x, "lower") + get_index(d, rep(0, length(attr(x, "lower"))),
      attr(x, "upper") - attr(x, "lower"))
    out[[d]]$data.name <- deparse_varname(attr(x, "name"), ind, ind)

    class(out[[d]]) <- "density.mcmcarray.univariate"
  }

  dim(out) <- dimen[!drop_dim]
  class(out) <- "density.mcmcarray"
  return(out)
}


#' @export
#' @rdname mcmcarray-object
biips_hist <- function(x, ...) UseMethod("biips_hist")

#' @export
#' @rdname mcmcarray-object
biips_hist.mcmcarray <- function(x, main = NULL, xlab = NULL, ...) {
  stopifnot(is.mcmcarray(x))
  out <- list()

  dimen <- dim(x)
  drop_dim <- names(dimen) %in% c("iteration", "chain")
  n_samples <- prod(dimen[drop_dim])
  len <- prod(dimen[!drop_dim])

  for (d in 1:len) {
    ind_vec <- seq(d, len * (n_samples - 1) + d, len)
    values <- x[ind_vec]

    ind <- attr(x, "lower") + get_index(d, rep(0, length(attr(x, "lower"))),
      attr(x, "upper") - attr(x, "lower"))
    xname <- deparse_varname(attr(x, "name"), ind, ind)

    main_ <- rec(main, d)
    if (is.null(main_)) {
      main_ <- paste("Histogram of", xname)
    }

    out[[d]] <- hist(values, main = main_, xlab = rec(xlab, d), ...)  # recycle arguments
    out[[d]]$xname <- xname
  }

  dim(out) <- dimen[!drop_dim]
  class(out) <- "histogram.mcmcarray"
  return(invisible(out))
}

#' @export
#' @rdname mcmcarray-object
biips_table.mcmcarray.list <- function(x, ...) {
  stopifnot(is.mcmcarray.list(x))
  out <- list()
  for (i in 1:length(x)) {
    name <- names(x)[i]
    if (!is.mcmcarray(x[[i]]) || name %in% c("log_marg_like_pen", "log_marg_like",
      "info"))
      next
    out[[name]] <- biips_table(x[[i]], ...)
  }
  class(out) <- "table.mcmcarray.list"
  return(out)
}

#' @export
#' @rdname mcmcarray-object
biips_density.mcmcarray.list <- function(x, bw = "nrd0", ...) {
  stopifnot(is.mcmcarray.list(x))
  out <- list()
  for (i in 1:length(x)) {
    name <- names(x)[i]
    if (!is.mcmcarray(x[[i]]) || name %in% c("log_marg_like_pen", "log_marg_like",
      "info"))
      next
    out[[name]] <- biips_density(x[[i]], bw = rec(bw, i), ...)  # recycle bw
  }
  class(out) <- "density.mcmcarray.list"
  return(out)
}

#' @export
#' @rdname mcmcarray-object
#' @param main,xlab plotting parameters with useful defaults.
biips_hist.mcmcarray.list <- function(x, main = NULL, xlab = NULL, ...) {
  stopifnot(is.mcmcarray.list(x))
  out <- list()
  for (i in 1:length(x)) {
    name <- names(x)[i]
    if (!is.mcmcarray(x[[i]]) || name %in% c("log_marg_like_pen", "log_marg_like",
      "info"))
      next
    out[[name]] <- biips_hist(x[[i]], main = rec(main, i), xlab = rec(xlab, i),
      ...)  # recycle arguments
  }
  class(out) <- "histogram.mcmcarray.list"
  return(invisible(out))
}


#' @export
#' @rdname mcmcarray-object
#' @return The method \code{summary} is an alias for \code{biips_summary}.
summary.mcmcarray <- function(object, ...) {
  return(biips_summary(object, ...))
}
#' @export
#' @rdname mcmcarray-object
summary.mcmcarray.list <- function(object, ...) {
  return(biips_summary(object, ...))
}

#' @importFrom stats density
#' @export
#' @rdname mcmcarray-object
#' @return The method \code{density} is an alias for \code{biips_density}.
#' @seealso \code{\link[stats]{density}}
density.mcmcarray <- function(x, ...) {
  return(biips_density(x, ...))
}
#' @export
#' @rdname mcmcarray-object
density.mcmcarray.list <- function(x, ...) {
  return(biips_density(x, ...))
}
#' @importFrom graphics hist
#' @export
#' @rdname mcmcarray-object
#' @return The method \code{hist} is an alias for \code{biips_hist}.
hist.mcmcarray <- function(x, ...) {
  return(biips_hist(x, ...))
}
#' @export
#' @rdname mcmcarray-object
hist.mcmcarray.list <- function(x, ...) {
  return(biips_hist(x, ...))
}
