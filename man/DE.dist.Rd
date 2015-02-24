\name{DE.dist}
\alias{DE.dist}
\title{Generates \code{iterations} randomly generated models then fits them to the
dataset}
\usage{
  DE.dist(data = NULL, paths, location = "",
    var.names = NULL, restrictions = NULL, arrows = 0.2,
    actual.model = NULL, actual.fit = NULL, iterations = 1,
    N = 1000, graphviz = FALSE, output.ram = TRUE,
    openmx = FALSE, check.identical.models = TRUE,
    limit = 100, upper.limit = TRUE, recursive = TRUE,
    fix.var = FALSE, fitIndex = c("RMSEA", "RMSR"),
    corr.exo = FALSE, corr.residuals = 0, ...)
}
\arguments{
  \item{data}{a correlation matrix containing the data to
  be fit}

  \item{paths}{the number of paths for each randomly
  generated model. Can either be a scalar, or a vector of
  length two, where the first element is the minimum number
  of paths and the second element is the maximum. See
  examples.}

  \item{location}{file location where the results are to be
  stored (as either a pdf or a RAM matrix)}

  \item{var.names}{the names of the variables}

  \item{restrictions}{a matrix containing restrictions on
  the DE procedure. See details section.}

  \item{arrows}{specifies what proportion of the paths
  should be bi-directional (correlational)}

  \item{actual.model}{a RAM matrix (containing columns
  from, to, arrows, values) of the proposed model. If
  specified, either the point estimate of RMSEA or its
  upper limit is used to identify confounding models.}

  \item{actual.fit}{the RMSEA of the actual model. Serves
  as a threshold for when a random model is outputted.}

  \item{iterations}{Number of models randomly generated.}

  \item{N}{the sample size that generated the correlation
  matrix.}

  \item{graphviz}{logical. When a model performs better
  than the specified threshold, should a graphic of the
  model be outputted? Requires prior installation of
  graphviz.}

  \item{output.ram}{logical. When a model performs better
  than the specified threshold, should a csv of the model
  be outputted?}

  \item{openmx}{logical. Should \code{OpenMx} be used to
  fit the model? If FALSE, the \code{sem} package is used
  instead. The default is FALSE since OpenMx cannot be
  installed from CRAN.}

  \item{check.identical.models}{logical. Should all models
  be unique? Increases computational intensity, but
  prevents duplicates in the DE procedure.}

  \item{limit}{The number of times the DE procedure should
  try to find a unique model before quitting early.
  Defaults at 100.}

  \item{recursive}{Should paths be limited to recursive
  paths?}

  \item{upper.limit}{logical. Should the upper limit of a
  90\% confidence interval for RMSEA be used to find
  confounding models?}

  \item{...}{Other arguments passed to the DE.graphviz
  function.}

  \item{fix.var}{Should the variances of the exogenous
  variables be fixed?}

  \item{fitIndex}{specifies which fit index to use for
  comparison. One of the following, either "RMSEA" or
  "RMSR."}

  \item{corr.exo}{Should the exogenous variables be
  automatically correlated?}

  \item{corr.residuals}{The proportion of arrows that
  should correlate residuals. Defaults to 0.}
}
\value{
  \item{DF}{The number of degrees of freedom for each
  randomly generated model.}

  \item{better}{An integer indicating the number of models
  that fit better}

  \item{fitIndex}{A string (either "RMSEA" or "RMSR")
  indicating which fit index was used.}
}
\description{
  This function randomly generates \code{iterations}
  models, where the user specifies the value of
  \code{iterations}. Each of these models are fit to the
  dataset provided and the fit is assessed. For those
  models that exceed a threshold, the model is outputted as
  a pdf graphic and/or csv RAM matrix.
}
\details{
  Model Conditioned Data Elasticity (DE) is a procedure
  proposed by Fife, Rodgers, and Mendoza (2014). Given a
  correlation matrix, the procedure randomly generates path
  analysis models. Those models that fit better than the
  one proposed are flagged (either outputted as a path
  diagram via graphviz, or outputted as a RAM matrix in csv
  form) for later inspection.

  The DE procedure also allows user restrictions. For
  example, suppose one variable is Age. Obviously, Age
  should probably not be endogenous, so the user can
  specify Age as an endogenous variable. That is done by
  creating a matrix where the columns correspond to "From",
  "To", and "Include." For example, to specify that A must
  cause B, one would insert in the first row of the matrix
  c("A", "B", "1"). To specify that nothing can cause a
  variable (i.e., to make a variable exogenous), one would
  leave the "From" column as "". For example, the Age
  example would have c("", "Age", "0").
}
\examples{
#### load the albanese dataset
data(albanese)

	#### set restrictions (make Age exogenous)
restrictions = matrix(c("", "Age", 0), nrow=1, byrow=TRUE)

	#### make a RAM matrix of the albanese model
alb = matrix(c("Age", "External", 1, .2,
				"Age", "Reflexive", 1, .2,
				"Age", "CPM", 1, .2,
				"CPM", "Reflexive", 1, .2,
				"CPM", "Mental", 1, .2,
				"Age", "Mental", 1,  .2), nrow=6, byrow=TRUE)
alb = data.frame(alb); names(alb) = c("From", "To", "Arrows", "Values")

	#### do DE
FF = DE.dist(data=albanese, paths=c(6,7), N=366, restrictions=restrictions,
		location="deleteme", var.names=names(albanese),
		actual.model=alb, iterations=20, arrows=0, openmx=TRUE,
		upper.limit=TRUE, iterations=50)
}
\author{
  Dustin Fife
}
\references{
  Fife, D.A., Rodgers, J.L., & Mendoza, J. L. (2013). Model
  conditioned data elasticity in path analysis: Assessing
  the "confoundability" of the data. Manuscript submitted
  for publication (under second revision).
}
\seealso{
  See Also \code{\link{DE}}, \code{\link{DE.fit}},
  \code{\link{print.DE.dist}}
}

