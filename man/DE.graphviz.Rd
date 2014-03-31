\name{DE.graphviz}
\alias{DE.graphviz}
\title{Outputs a RAM model as a graphic}
\usage{
  DE.graphviz(model, file = "deleteme",
    dotFilename = "graphic", rank.direction = "LR",
    max.rank = NULL, min.rank = NULL, same.rank = NULL,
    delete.dot = FALSE, silent = TRUE, label = "",
    node.labels = TRUE)
}
\arguments{
  \item{model}{The RAM model to be imaged.}

  \item{file}{location where the pdf file is to be stored}

  \item{dotFilename}{not used.}

  \item{rank.direction}{a string of the form, either "LR",
  or "TB", indicating Left to right or top to bottom}

  \item{max.rank}{a string of the form "V1, V2, V3", where
  variables 1-3 are either at the most-right of the page or
  top-most (depending on whether rank.direction = LR or
  TB)}

  \item{min.rank}{same as above but will be in most-left or
  bottom-most}

  \item{same.rank}{what variables should be lined up. Uses
  same format as max.rank and min.rank.}

  \item{delete.dot}{logical. Indicates whether the dot
  files should be deleted (and thus only return pdf files)}

  \item{silent}{Should R output information from system
  commands?}

  \item{label}{The label attached to the graphviz diagram.
  It will be displayed at the bottom of the image.}

  \item{node.labels}{Should the nodes be labeled?}
}
\description{
  Outputs a pdf graphic of the model using graphviz
}
\examples{
# generate random model
rand.model = DE(LETTERS[1:5], paths=6)
# output to graphic
DE.graphviz(rand.model, min.rank="A, B, C")
# now go look at the root directory for the pdf graphic
}
\author{
  Dustin Fife
}

