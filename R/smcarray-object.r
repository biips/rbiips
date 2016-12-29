#' @title Objects for representing SMC output.
#'
#' @name smcarray-object
#' @aliases smcarray.fsb-object smcarray.fsb.list-object smcarray smcarray.fsb
#'   smcarray.fsb.list
#'
#' @description A \code{smcarray} object is used by the
#' \code{\link{biips_smc_samples}} function to represent SMC output or particles
#' of a given variable.
#'
#' A \code{smcarray.fsb} object is a named list of \code{smcarray} objects with
#' different types of monitoring for the same variable. Members in this list
#' have names \code{f} (filtering), \code{s} (smoothing) or \code{b} (backward
#' smoothing).
#'
#' A \code{smcarray.fsb.list} object is a named list of \code{smcarray.fsb}
#' objects for different monitored variables. It might also contain a member
#' named \code{log_marg_like} with an estimate of the log marginal likelihood.
#'
#' The methods apply identically to \code{smcarray}, \code{smcarray.fsb} or
#' \code{smcarray.fsb.list} objects and return a named list with the same named
#' members as the input object.
#'
#' @details Assuming \code{dim} is the dimension of the monitored variable, a
#'   \code{smcarray} object is a list with the members:
#'   \describe{
#'     \item{values}{array of dimension \code{c(dim, n_part)} with the values of
#'       the particles.}
#'     \item{weights}{array of dimension \code{c(dim, n_part)} with the weights
#'       of the particles.}
#'     \item{ess}{array of dimension \code{dim} with Effective Sample Sizes of
#'       the particles set.}
#'     \item{discrete}{array of dimension \code{dim} with logicals indicating
#'       discreteness of each component.}
#'     \item{iterations}{array of dimension \code{dim} with sampling iterations
#'       of each component.}
#'     \item{conditionals}{lists of the contitioning variables (observations). Its
#'       value is:
#'       \itemize{
#'         \item for filtering: a list of dimension \code{dim}. each member is a
#'           character vector with the respective conditioning variables of
#'           the node array component.
#'         \item for smoothing/backward_smoothing: a character vector, the same
#'           for all the components of the node array.
#'       }}
#'     \item{name}{string with the name of the variable.}
#'     \item{lower}{vector with the lower bounds of the variable.}
#'     \item{upper}{vector with the upper bounds of the variable.}
#'     \item{type}{string with the type of monitor (\code{'filtering'},
#'       \code{'smoothing'} or \code{'backward_smoothing'}).}
#'   }
#'
#'   For instance, if \code{out_smc} is a \code{smcarray.fsb.list} object, one can
#'   access the values of the smoothing particles for the variable \code{'x'}
#'   with: \code{out_smc$x$s$values}.
#'
#' @param object,x a \code{smcarray}, \code{smcarray.fsb} or
#'   \code{smcarray.fsb.list} object.
#' @param ... additional arguments to be passed to the default methods. See
#'   \code{\link[stats]{density}}, \code{\link[base]{table}}
#'
#' @return The methods apply identically to \code{smcarray}, \code{smcarray.fsb} or
#'   \code{smcarray.fsb.list} objects and return a named list with the same
#'   named members as the input object.
#'
#' @examples
#' modelfile <- system.file('extdata', 'hmm.bug', package = 'rbiips')
#' stopifnot(nchar(modelfile) > 0)
#' cat(readLines(modelfile), sep = '\n')
#'
#' data <- list(tmax = 10, p = c(.5, .5), logtau_true = log(1), logtau = log(1))
#' model <- biips_model(modelfile, data, sample_data = TRUE)
#'
#' n_part <- 100
#' out_smc <- biips_smc_samples(model, c('x', 'c[2:10]'), n_part, type = 'fs',
#'                              rs_thres = 0.5, rs_type = 'stratified')
#'
#' #' Manipulate `smcarray.fsb.list` object
#' is.smcarray.fsb.list(out_smc)
#' names(out_smc)
#' out_smc
#' biips_diagnosis(out_smc)
#' biips_summary(out_smc)
#'
#' #' Manipulate `smcarray.fsb` object
#' is.smcarray.fsb(out_smc$x)
#' names(out_smc$x)
#' out_smc$x
#' biips_diagnosis(out_smc$x)
#' summ_smc_x <- biips_summary(out_smc$x, order = 2, probs = c(.025, .975))
#' summ_smc_x
#' dens_smc_x <- biips_density(out_smc$x, bw = 'nrd0', adjust = 1, n = 100)
#' par(mfrow = c(2, 2))
#' plot(dens_smc_x)
#'
#' is.smcarray.fsb(out_smc[['c[2:10]']])
#' names(out_smc[['c[2:10]']])
#' out_smc[['c[2:10]']]
#' biips_diagnosis(out_smc[['c[2:10]']])
#' summ_smc_c <- biips_summary(out_smc[['c[2:10]']])
#' summ_smc_c
#' table_smc_c <- biips_table(out_smc[['c[2:10]']])
#' par(mfrow = c(2, 2))
#' plot(table_smc_c)
#'
#' #' Manipulate `smcarray` object
#' is.smcarray(out_smc$x$f)
#' names(out_smc$x$f)
#' out_smc$x$f
#' out_smc$x$s
#' biips_diagnosis(out_smc$x$f)
#' biips_diagnosis(out_smc$x$s)
#' biips_summary(out_smc$x$f)
#' biips_summary(out_smc$x$s)
#' par(mfrow = c(2, 2))
#' plot(biips_density(out_smc$x$f))
#' par(mfrow = c(2, 2))
#' plot(biips_density(out_smc$x$s))
#'
#' par(mfrow = c(2, 2))
#' plot(biips_table(out_smc[['c[2:10]']]$f))
#' par(mfrow = c(2, 2))
#' plot(biips_table(out_smc[['c[2:10]']]$s))
NULL


#' @export
#' @rdname smcarray-object
#' @return The function \code{is.smcarray} returns \code{TRUE} if the object is of class \code{smcarray}.
is.smcarray <- function(object) {
  return(class(object) == "smcarray")
}
#' @export
#' @rdname smcarray-object
#' @return The function \code{is.smcarray.fsb} returns \code{TRUE} if the object
#'   is of class \code{smcarray.fsb}.
is.smcarray.fsb <- function(object) {
  return(class(object) == "smcarray.fsb")
}
#' @export
#' @rdname smcarray-object
#' @return The function \code{is.smcarray.fsb.list} returns \code{TRUE} if the
#'   object is of class \code{smcarray.fsb.list}.
is.smcarray.fsb.list <- function(object) {
  return(class(object) == "smcarray.fsb.list")
}

#' @export
print.smcarray <- function(x, ...) {
  print(summary(x), ...)
}
#' @export
print.smcarray.fsb <- function(x, ...) {
  print(summary(x), ...)
}
#' @export
print.smcarray.fsb.list <- function(x, ...) {
  print(summary(x), ...)
  if ("log_marg_like" %in% names(x))
    cat("Log-marginal likelihood: ", x$log_marg_like, "\n")
}

#' @export
#' @rdname smcarray-object
biips_diagnosis <- function(object, ...) UseMethod("biips_diagnosis")

#' @export
#' @rdname smcarray-object
#' @param ess_thres  integer. Threshold on the Effective Sample Size (ESS). If
#'   all the ESS components are over \code{ess_thres}, the diagnostic is
#'   \code{'GOOD'}, otherwise it is \code{'BAD'}. (default=30).
#' @param quiet  logical. Disable message display. (default=\code{FALSE}).
#'
#' @return The method \code{biips_diagnosis} prints diagnosis of the SMC output
#'   and returns the minimum ESS value.
biips_diagnosis.smcarray <- function(object, ess_thres = 30, quiet = FALSE, ...) {
  stopifnot(is.smcarray(object))
  stopifnot(is.numeric(ess_thres), length(ess_thres) == 1, is.finite(ess_thres),
    ess_thres >= 0)
  stopifnot(is.logical(quiet), length(quiet) == 1)

  ess_min <- min(object$ess)
  out <- list(ess_min = ess_min, valid = ess_min > ess_thres)
  class(out) <- "diagnosis.smcarray"

  if (!quiet) {
    varname <- deparse_varname(object$name, object$lower, object$upper)
    cat("* Diagnosis of variable:", varname, "\n")
    switch(object$type, filtering = {
      cat("  Filtering: ")
    }, smoothing = {
      cat("  Smoothing: ")
    }, backward_smoothing = {
      cat("  Backward smoothing: ")
    })
    print(out)
  }

  return(invisible(out))
}

#' @export
#' @rdname smcarray-object
#' @param type string containing the characters \code{'f'} (fitering),
#'   \code{'s'} (smoothing) and/or \code{'b'} (backward smoothing).
#'   Select the corresponding members of the input to be analysed.
#'   (default = \code{'fsb'}).
biips_diagnosis.smcarray.fsb <- function(object, type = "fsb", quiet = FALSE, ...) {
  stopifnot(is.smcarray.fsb(object))
  type <- check_type(type)

  if (!quiet) {
    varname <- deparse_varname(object[[1]]$name, object[[1]]$lower, object[[1]]$upper)
    cat("* Diagnosis of variable:", varname, "\n")
  }

  out <- list()
  for (n in names(object)) {
    if (!(n %in% type))
      next
    out[[n]] <- biips_diagnosis(object[[n]], quiet = TRUE, ...)
    if (!quiet) {
      switch(n, f = {
        cat("  Filtering: ")
      }, s = {
        cat("  Smoothing: ")
      }, b = {
        cat("  Backward smoothing: ")
      })
      print(out[[n]])
    }
  }

  return(invisible(out))
}

#' @export
#' @rdname smcarray-object
biips_diagnosis.smcarray.fsb.list <- function(object, type = "fsb", quiet = FALSE,
  ...) {
  stopifnot(is.smcarray.fsb.list(object))

  out <- list()
  for (n in names(object)) {
    if (!is.smcarray.fsb(object[[n]]))
      next
    out <- biips_diagnosis(object[[n]], type = type, quiet = quiet, ...)
  }

  return(invisible(out))
}

#' @export
print.diagnosis.smcarray <- function(x, ...) {
  if (x$valid) {
    cat("GOOD\n")
  } else {
    cat("POOR\n")
    cat("    The minimum effective sample size is too low:", x$ess_min, "\n")
    cat("    Estimates may be poor for some variables.\n")
    cat("    You should increase the number of particles\n.")
  }
  invisible()
}


#' @export
#' @rdname smcarray-object
biips_summary <- function(object, ...) UseMethod("biips_summary")

#' @rdname smcarray-object
#' @export
#' @param probs    vector of reals. probability levels in ]0,1[ for quantiles.
#'   (default = \code{c()})
#' @param order    integer. Moment statistics of order below or equal to
#'   \code{order} are returned. (default = 0 if all the components are discrete
#'   variables and 1 otherwise)
#' @param mode     logical. Activate computation of the mode, i.e. the most
#'   frequent value among the particles. (default = \code{TRUE} if all the components
#'   are discrete variables and \code{FALSE} otherwise)
#' @return The method \code{biips_summary} returns univariate marginal
#'   statistics. The output innermost members are objects of class
#'   \code{summary.smcarray}. Assuming \code{dim} is the dimension of the
#'   variable, the \code{summary.smcarray} object is a list with the following
#'   members:
#'   \item{mean}{array of size \code{dim}. The mean if \code{order>=1}.}
#'   \item{var}{array of size \code{dim}. The variance, if \code{order>=2}.}
#'   \item{skew}{array of size \code{dim}. The skewness, if \code{order>=3}.}
#'   \item{kurt}{array of size \code{dim}. The kurtosis, if \code{order>=4}.}
#'   \item{probs}{vector of quantile probabilities.}
#'   \item{quant}{list of arrays of size \code{dim} for each probability level
#'       in \code{probs}. The quantile values, if \code{probs} is not empty.}
#'   \item{mode}{array of size \code{dim}. The most frequent values for
#'       discrete components.}
biips_summary.smcarray <- function(object, probs = c(), order = ifelse(mode, 0, 1),
  mode = all(object$discrete), ...) {
  stopifnot(is.smcarray(object))

  ### TODO check arguments
  if (length(probs) > 0)
    stopifnot(is.numeric(probs), probs > 0, probs < 1)

  drop_dims <- names(dim(object$values)) %in% c("particle")
  n_part <- dim(object$values)["particle"]
  dim_array <- dim(object$values)[!drop_dims]
  len <- prod(dim_array)
  m <- which(!drop_dims)
  summ <- list()

  ### moment based statistics
  if (order >= 1)
    summ$mean <- apply(object$values * object$weights, m, FUN = sum)

  if (order >= 2)
    summ$var <- apply(mapply("-", object$values, summ$mean)^2 * object$weights,
      m, FUN = sum)

  if (order >= 3) {
    mom2 <- apply(object$values^2 * object$weights, m, FUN = sum)
    mom3 <- apply(object$values^3 * object$weights, m, FUN = sum)
    summ$skew <- (mom3 - 3 * mom2 * summ$mean + 2 * summ$mean^3)/summ$var^(3/2)
  }

  if (order >= 4) {
    mom4 <- apply(object$values^4 * object$weights, m, FUN = sum)
    summ$kurt <- (mom4 - 4 * mom3 * summ$mean + 6 * mom2 * summ$mean^2 - 3 *
      summ$mean^4)/summ$var^2 - 3
  }

  ### quantiles
  if (length(probs) > 0) {
    summ$probs <- probs
    summ$quant <- list()

    for (d in 1:len) {
      indvec <- seq(d, len * (n_part - 1) + d, len)
      stat_d <- rbiips("wtd_quantile", object$values[indvec], n_part * object$weights[indvec],
        probs)
      quant_names <- names(stat_d)
      if (d == 1) {
        for (n in quant_names) {
          if (is.null(summ$quant[[n]]))
          summ$quant[[n]] <- array(dim = dim_array)
        }
      }
      for (n in quant_names) summ$quant[[n]][d] <- stat_d[[n]]
    }
  }

  ### mode
  if (mode) {
    summ$mode <- array(dim = dim_array)

    for (d in 1:len) {
      indvec <- seq(d, len * (n_part - 1) + d, len)
      stat_d <- rbiips("wtd_mode", object$values[indvec], n_part * object$weights[indvec])
      summ$mode[d] <- stat_d
    }
  }

  summ$drop_dims <- dim(object$values)[drop_dims]

  class(summ) <- "summary.smcarray"

  return(summ)
}

#' @export
#' @rdname smcarray-object
biips_summary.smcarray.fsb <- function(object, ...) {
  stopifnot(is.smcarray.fsb(object))
  out <- list()
  for (n in names(object)) out[[n]] <- biips_summary(object[[n]], ...)

  class(out) <- "summary.smcarray.fsb"

  return(out)
}

#' @export
#' @rdname smcarray-object
biips_summary.smcarray.fsb.list <- function(object, ...) {
  stopifnot(is.smcarray.fsb.list(object))
  out <- list()

  for (n in names(object)) {
    if (!is.smcarray.fsb(object[[n]]))
      next
    out[[n]] <- biips_summary(object[[n]], ...)
  }

  class(out) <- "summary.smcarray.fsb.list"

  return(out)
}


#' @export
print.summary.smcarray <- function(x, ...) {
  cat("smcarray:\n")
  print(x[!(names(x) %in% c("drop_dims"))], ...)
  if (length(x$drop_dims) > 0) {
    cat("Marginalizing over:", paste(paste(names(x$drop_dims), "(", x$drop_dims,
      ")", sep = ""), collapse = ","), "\n")
  }
  invisible()
}

#' @export
print.summary.smcarray.fsb <- function(x, prefix = NULL, ...) {
  for (n in names(x)) {
    if (!is.null(prefix))
      cat(prefix, " ", sep = "")
    switch(n, f = cat("filtering "), s = cat("smoothing "), f = cat("backward smoothing "))
    print(x[[n]], ...)
    cat("\n")
  }
  invisible()
}

#' @export
print.summary.smcarray.fsb.list <- function(x, ...) {
  for (n in names(x)) {
    print(x[[n]], prefix = n, ...)
    cat("\n")
  }
  invisible()
}


#' @export
#' @rdname smcarray-object
biips_table <- function(x, ...) UseMethod("biips_table")

#' @export
#' @rdname smcarray-object
#' @return The method \code{biips_table} returns univariate marginal frequency
#'   tables or probability mass estimates of discrete variables. The output
#'   innermost members are objects of class \code{table.smcarray}.
biips_table.smcarray <- function(x, ...) {
  stopifnot(is.smcarray(x))
  out <- list()

  dimen <- dim(x$values)
  drop_dim <- names(dimen) %in% c("particle")

  n_part <- dimen["particle"]
  if (is.na(n_part))
    n_part <- 1

  len <- prod(dimen[!drop_dim])

  for (d in 1:len) {
    ind_vec <- seq(d, len * (n_part - 1) + d, len)
    values <- x$values[ind_vec]
    weights <- x$weights[ind_vec]

    out[[d]] <- rbiips("wtd_table", values, weights)
    dn <- list(names(out[[d]]))
    ind <- x$lower + get_index(d, rep(0, length(x$lower)), x$upper - x$lower)
    names(dn) <- deparse_varname(x$name, ind, ind)

    out[[d]] <- array(rbiips("wtd_table", values, weights), dimnames = dn)
    out[[d]] <- as.table(out[[d]])
  }

  dim(out) <- dimen[!drop_dim]
  class(out) <- "table.smcarray"
  return(out)
}



#' @export
#' @rdname smcarray-object
biips_density <- function(x, ...) UseMethod("biips_density")

#' @export
#' @rdname smcarray-object
#' @param bw   either a real with the smoothing bandwidth to be used or a string
#'   giving a rule to choose the bandwidth. See \code{\link[stats]{bw.nrd}}.
#'   (default=\code{'nrd0'})
#' @return The method \code{biips_density} returns univariate marginal kernel
#'   density estimates. The output innermost members are objects of class
#'   \code{density.smcarray}.
biips_density.smcarray <- function(x, bw = "nrd0", ...) {
  stopifnot(is.smcarray(x))
  out <- list()

  dimen <- dim(x$values)
  drop_dim <- names(dimen) %in% c("particle")

  n_part <- dimen["particle"]
  if (is.na(n_part))
    n_part <- 1

  len <- prod(dimen[!drop_dim])

  for (d in 1:len) {
    ind_vec <- seq(d, len * (n_part - 1) + d, len)
    values <- x$values[ind_vec]
    weights <- x$weights[ind_vec]

    out[[d]] <- density(values, weights = weights, bw = rec(bw, d), ...)  # recycle bw

    ind <- x$lower + get_index(d, rep(0, length(x$lower)), x$upper - x$lower)
    out[[d]]$data.name <- deparse_varname(x$name, ind, ind)
  }

  dim(out) <- dimen[!drop_dim]
  class(out) <- "density.smcarray"
  return(out)
}

#' @export
#' @rdname smcarray-object
biips_table.smcarray.fsb <- function(x, ...) {
  stopifnot(is.smcarray.fsb(x))
  out <- list()

  for (fsb in names(x)) {
    if (!is.smcarray(x[[fsb]]))
      next
    dimen <- dim(x[[fsb]]$values)
    tab <- biips_table(x[[fsb]], ...)
    for (d in 1:length(tab)) {
      if (length(out) < d)
        out[[d]] <- list()
      out[[d]][[fsb]] <- tab[[d]]
    }
  }

  for (d in 1:length(tab)) class(out[[d]]) <- "table.smcarray.fsb.univariate"

  drop_dim <- names(dimen) %in% c("particle")
  dim(out) <- dimen[!drop_dim]

  class(out) <- "table.smcarray.fsb"
  return(out)
}

#' @export
#' @rdname smcarray-object
#' @param adjust  scale factor for the bandwidth. the bandwidth used is actually
#'   \code{adjust*bw}. This makes it easy to specify values like \sQuote{half
#'   the default} bandwidth. (default = 1)
biips_density.smcarray.fsb <- function(x, bw = "nrd0", adjust = 1, ...) {
  stopifnot(is.smcarray.fsb(x))
  out <- list()
  bw_s <- bw

  # first treat filtering and backward smoothing
  if (!is.null(x$f)) {
    dimen <- dim(x$f$values)
    dens <- biips_density(x$f, bw, adjust = adjust, ...)
    bw_s <- list()
    for (d in 1:length(dens)) {
      if (length(out) < d)
        out[[d]] <- list()
      out[[d]]$f <- dens[[d]]
      bw_s[[d]] <- dens[[d]]$bw
    }
  }
  if (!is.null(x$b)) {
    dimen <- dim(x$b$values)
    dens <- biips_density(x$b, bw, adjust = adjust, ...)
    bw_s <- list()
    for (d in 1:length(dens)) {
      if (length(out) < d)
        out[[d]] <- list()
      out[[d]]$b <- dens[[d]]
      bw_s[[d]] <- dens[[d]]$bw
    }
  }

  # then treat smoothing (applying previous bandwidth if any)
  if (!is.null(x$s)) {
    dimen <- dim(x$s$values)
    if (is.list(bw_s))
      adjust <- 1
    dens <- biips_density(x$s, bw_s, adjust = adjust, ...)
    for (d in 1:length(dens)) {
      if (length(out) < d)
        out[[d]] <- list()
      out[[d]]$s <- dens[[d]]
    }
  }

  for (d in 1:length(dens)) class(out[[d]]) <- "density.smcarray.fsb.univariate"

  drop_dim <- names(dimen) %in% c("particle")
  dim(out) <- dimen[!drop_dim]

  class(out) <- "density.smcarray.fsb"
  return(out)
}


#' @export
#' @rdname smcarray-object
biips_table.smcarray.fsb.list <- function(x, ...) {
  stopifnot(is.smcarray.fsb.list(x))
  out <- list()
  for (i in 1:length(x)) {
    if (!is.smcarray.fsb(x[[i]]))
      next
    name <- names(x)[i]
    out[[name]] <- biips_table(x[[i]], ...)
  }

  class(out) <- "table.smcarray.fsb.list"
  return(out)
}


#' @export
#' @rdname smcarray-object
biips_density.smcarray.fsb.list <- function(x, bw = "nrd0", ...) {
  stopifnot(is.smcarray.fsb.list(x))
  out <- list()
  for (i in 1:length(x)) {
    if (!is.smcarray.fsb(x[[i]]))
      next
    name <- names(x)[i]
    out[[name]] <- biips_density(x[[i]], bw = rec(bw, i), ...)  # recycle bw
  }

  class(out) <- "density.smcarray.fsb.list"

  return(out)
}


#' @export
#' @rdname smcarray-object
#' @return The method \code{summary} is an alias for \code{biips_summary}.
summary.smcarray <- function(object, ...) {
  return(biips_summary(object, ...))
}
#' @export
#' @rdname smcarray-object
summary.smcarray.fsb <- function(object, ...) {
  return(biips_summary(object, ...))
}
#' @export
#' @rdname smcarray-object
summary.smcarray.fsb.list <- function(object, ...) {
  return(biips_summary(object, ...))
}

#' @importFrom stats density
#' @rdname smcarray-object
#' @export
#' @return The method \code{density} is an alias for \code{biips_density}.
#' @seealso \code{\link[stats]{density}}
density.smcarray <- function(x, ...) {
  return(biips_density(x, ...))
}
#' @export
#' @rdname smcarray-object
density.smcarray.fsb <- function(x, ...) {
  return(biips_density(x, ...))
}
#' @export
#' @rdname smcarray-object
density.smcarray.fsb.list <- function(x, ...) {
  return(biips_density(x, ...))
}
