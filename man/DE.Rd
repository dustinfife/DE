\name{DE}
\alias{DE}
\title{Generates a random model using RAM specification.}
\usage{
  DE(variable.names, paths, restrictions = NULL,
    prop.arrows = 0.2, allow.orphaned = FALSE,
    allow.bidir = FALSE, corr.exogenous = FALSE,
    corr.residuals = 0)
}
\arguments{
  \item{variable.names}{A vector of variable names.}

  \item{paths}{The number of paths to be randomly
  generated. Can be a single value or a vector of two
  integers specifying a range.}

  \item{restrictions}{What kind of restrictions are set.
  (See details).}

  \item{prop.arrows}{What proportion of exogenous variables
  should be correlated. Defaults to .2.}

  \item{allow.orphaned}{Should orphaned variables be
  allowed when random models are generated?}

  \item{allow.bidir}{Should bidirectional arrows be
  allowed? (Note: this is not the same as a correlation).}

  \item{corr.exogenous}{Should exogenous variables be
  correlated?}

  \item{corr.residuals}{a value between 0 and 1 indicating
  the proportion of residuals the user allows to be
  correlated}
}
\value{
  Returns a RAM matrix.
}
\description{
  This function generates a single random path model.
}
\details{
  Making restrictions is simple. For example, suppose one
  variable is Age. Obviously, Age should probably not be
  endogenous, so the user can specify Age as an endogenous
  variable. That is done by creating a matrix where the
  columns correspond to "From", "To", and "Include." For
  example, to specify that A must cause B, one would insert
  in the first row of the matrix c("A", "B", "1"). To
  specify that nothing can cause a variable (i.e., to make
  a variable exogenous), one would leave the "From" column
  as "". For example, the Age example would have c("",
  "Age", "0").

  Allowing any variable to correlate with an endogenous
  variable is equivalent to correlating with the residuals
  of that endogenous variable. When the user specifies a
  non-zero value (k) for corr.residuals, the algorithm
  randomly selects k*(number of paths) of the paths to be
  double-headed, thereby permitting correlated residuals.
}
\examples{
rest = matrix(c("A", "B", "1", 1,
						"", "A", "0", 1), nrow=2, byrow=TRUE)
DE(variable.names=LETTERS[1:6], paths=c(6,7), restrictions = NULL, prop.arrows = 0.2, allow.orphaned=FALSE, allow.bidir=FALSE, allow.cov.endogenous=FALSE)
}
\author{
  Dustin Fife
}
\references{
  Fife, D.A., Rodgers, J.L., & Mendoza, J. L. (2013). Model
  conditioned data elasticity in path analysis: Assessing
  the "confoundability"
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}

