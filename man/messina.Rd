\name{messina}
\alias{messina}
\title{Identify features that can optimally separate two classes of samples, subject
to minimum performance requirements.}
\usage{
  messina(x, y, min_sens, min_spec, f_train = 0.9,
    n_boot = 50, seed = NULL)
}
\arguments{
  \item{x}{feature expression values, either supplied as an
  ExpressionSet, or as an object that can be converted to a
  matrix by as.matrix.  In the latter case, features should
  be in rows and samples in columns, with feature names
  taken from the rows of the object.}

  \item{y}{a binary vector (TRUE/FALSE or 1/0) of class
  membership information for each sample in x.}

  \item{min_sens}{the minimum acceptable sensitivity that a
  classifier separating the two groups of y must achieve.}

  \item{min_spec}{the minimum acceptable specificity that a
  classifier separating the two groups of y must achieve.}

  \item{f_train}{the fraction of samples to be used in the
  training splits of the bootstrap rounds.}

  \item{n_boot}{the number of bootstrap rounds to use.}

  \item{seed}{an optional random seed for the analysis.}
}
\value{
  an object of class "MessinaResult" containing the results
  of the analysis.
}
\description{
  TODO
}
\author{
  Mark Pinese \email{m.pinese@garvan.org.au}
}
\references{
  TODO
}
\seealso{
  \code{\link{MessinaResult-class}}

  \code{\link{ExpressionSet}}
}

