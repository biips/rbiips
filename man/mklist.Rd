% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{mklist}
\alias{mklist}
\title{Make a list using names}
\usage{
mklist(names, env = parent.frame())
}
\arguments{
\item{names}{character strings of names of objects}

\item{env}{the environment to look for objects with names.}
}
\description{
Make a list using names
}
\note{
we use inherits = TRUE when calling mget
}
\author{
Adapted from rstan.
}
\keyword{internal}
