\name{DE.crit.path}
\alias{DE.crit.path}
\title{Find critical paths in a Path Analysis model.}
\usage{
DE.crit.path(data, variable.names = names(data), paths, prop.arrows = 0.2,
  iterations = 50, restrictions = NULL, cutpoint)
}
\arguments{
  \item{data}{The correlation matrix to be fit.}

  \item{variable.names}{A list of variables names if
  they're not specified in the correlation matrix.}

  \item{paths}{How many paths should be generated? Can be
  specified as a single value or as a range. See examples.}

  \item{prop.arrows}{What proportion of arrows should be
  correlational?}

  \item{iterations}{Number of random models to be
  generated.}

  \item{restrictions}{A matrix containing restrictions on
  the DE procedure. See details section.}

  \item{cutpoint}{An RMSEA value that differentiates
  between the modes of the two distributions.}
}
\value{
returns a matrix the specifies which paths are contained in
the "highMode" versus "lowMode" datasets
}
\description{
This function randomly generates models and determines
which paths, when removed, tend to reduce the fit of the
model substantially.
}
\details{
Fife, Rodgers, and Mendoza (2013) noted that DE
distributions are frequently bimodal and commented that a
possible reason for this is that the poor fitting
distribution has one or more ``critical paths'' that are
missing. This function first generates \code{iterations}
random models, then uses the value for \code{cutpoint} to
identify which paths exist in the right versus left side of
the cutpoint.
}
\examples{
restrictions = matrix(c("", "Age", 0), nrow=1, byrow=TRUE)
data(albanese)

crit = DE.crit.path(data=albanese, paths=c(5,7), prop.arrows=0, iterations = 20, restrictions=restrictions, cutpoint=.4)
crit
}
\author{
Dustin Fife
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

