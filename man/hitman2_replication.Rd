% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hitman2_replication.R
\name{hitman2_replication}
\alias{hitman2_replication}
\title{High-throughput approach for testing directional replication that uses filtering to improve adjusted p-values}
\usage{
hitman2_replication(
  tab,
  cols = 1:4,
  reorder.rows = FALSE,
  p.adj.rate = c("FDR", "FWER"),
  prefix = NA
)
}
\arguments{
\item{tab}{Matrix-like object with statistical and p-value columns. Only the signs of the statistics columns are used.
\code{tab} should have non-duplicated row names and should not have missing values.}

\item{cols}{Vector of column indices or names in the order of \code{c(stat1, p1, stat2, p2)}.}

\item{reorder.rows}{Logical, should rows be reordered by p-value?}

\item{p.adj.rate}{Either "FDR" for false discovery rate or "FWER" for family-wise error rate, the rate controlled by the Bonferroni procedure.}

\item{prefix}{Character string of length one with prefix of returned columns, e.g. if \code{prefix="repl"}, returned columns might
be \code{c("repl.chisq", "repl.p", "repl.FDR")}.
Prefix is not added if it is \code{NA}.}
}
\value{
Data frame whose rows correspond to the rows of \code{tab} with the same row names and whose columns are
\describe{
\item{chisq}{Chi-square for replication on 1 degreee of freedom.}
\item{p}{P-value for replication}
\item{FDR or FWER}{FDR or FWER for replication}
}
}
\description{
High-throughput approach for testing replication of two base studies, where each study applies a two-sided test whose
signed statistic and p-value are supplied via \code{tab}, and its desired to find rows where there is replication in a common direction.
Hitman2 improves on the adjusted p-values using filtering.
}
\details{
Larger chi-square values are more significant.
}
