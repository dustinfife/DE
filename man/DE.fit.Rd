\name{DE.fit}
\alias{DE.fit}
\title{Fit a RAM matrix to the data using path analysis.}
\usage{
  DE.fit(returned.model, dataset, N = 1000, openmx = FALSE,
    fix.variances = FALSE, calc.se = FALSE)
}
\arguments{
  \item{returned.model}{The model of interest, in RAM
  format. (A column of "From" variables, "To" variables,
  the number of arrows, and the starting value.)}

  \item{dataset}{A correlation matrix.}

  \item{N}{The sample sized used to estimate the
  correlation matrix.}

  \item{openmx}{logical. Should \code{OpenMx} be used to
  fit the model? If FALSE, the \code{sem} package is used
  instead. The default is FALSE since the \code{sem}
  package is faster.}

  \item{fix.variances}{logical. Should variances be fixed
  to one?}

  \item{calc.se}{Should the standard errors be calculated?
  Defaults to no to speed up the DE algorithm.}
}
\value{
  Depending on which program is used to fit the data,
  either an OpenMx or an sem output.
}
\description{
  Given a single RAM matrix, this function fits the model
  and returns the output from either \code{OpenMx} or
  \code{sem}.
}
\details{
  %% ~~ If necessary, more details than the description
  above ~~
}
\note{
  %% ~~further notes~~
}
\examples{
data(albanese)

	#### make restrictions
restrictions = matrix(c("", "Age", 0), nrow=1, byrow=TRUE)

	#### make a RAM matrix of the albanese model
alb = matrix(c("Age", "External", 1, .2,
				"Age", "Reflexive", 1, .2,
				"Age", "CPM", 1, .2,
				"CPM", "Reflexive", 1, .2,
				"CPM", "Mental", 1, .2,
				"Age", "Mental", 1,  .2), nrow=6, byrow=TRUE)
alb = data.frame(alb); names(alb) = c("From", "To", "Arrows", "Values")

model = DE.fit(returned.model=alb, dataset=albanese, N = 366, openmx=TRUE, fix.variances=TRUE)
summary(model)
}
\author{
  Dustin Fife
}
\references{
  %% ~put references to the literature/web site here ~
}
\seealso{
  \code{\link{DE}}, \code{\link{DE.dist}}
}
\keyword{~kwd1}
\keyword{~kwd2}

