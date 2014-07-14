\name{is.recur}
\alias{is.recur}
\title{Checks whether a model is recursive (i.e., limits to a DAG)}
\usage{
is.recur(r.model, var.names)
}
\arguments{
  \item{r.model}{The model generated with the function
  \code{\link{DE}}}

  \item{var.names}{The names of the variables}
}
\value{
Logical (either TRUE or FALSE)
}
\description{
Checks whether a model is recursive (i.e., limits to a DAG)
}
\seealso{
\code{\link{is.dag}}
}

