\name{ram.2.cor}
\alias{ram.2.cor}
\title{Generate a correlation matrix from a Reticular Action Model (RAM)}
\usage{
  ram.2.cor(RAM)
}
\arguments{
  \item{RAM}{a Reticular Action model with columns "From",
  "To", "Arrows", and "Values"}
}
\value{
  A correlation matrix generated from a RAM model
}
\description{
  Generate a correlation matrix from a RAM
}
\details{
  Although generating an implied covariance matrix from a
  RAM model is easy, doing so for a correlation is not as
  trivial since the diagonals are constrained to one. The
  symmetric matrix (S) must be formulated in such as way
  that the sum of the variances equals one. Although there
  may exist a closed-form solution, I know of none. So,
  instead, I simply used the optim function to do it. The
  user specifies the path coefficients in the form of a RAM
  matrix and the algorithm estimates what the residual
  correlations ought to be.
}
\examples{
## generate a simple exercise model
RAM = data.frame(matrix(c(
	"Age", "Energy", 1, .1,
	"Age", "Self-Efficacy", 2, .1,
	"Energy", "Exercise", 1, .6,
	"Eating", "Energy", 1, .4,
	"Eating", "Self-Efficacy", 1, .3,
	"Age", "Eating", 2, .15,
	"Self-Efficacy", "Exercise", 1, .5), ncol=4, byrow=T))
names(RAM) = c("From", "To", "Arrows", "Values")
ram.2.cor(RAM)
}
\author{
  Dustin Fife
}
\references{
  An openmx forum about this:
  http://openmx.psyc.virginia.edu/thread/3866
}

