\name{wine}
\alias{wine}
\docType{data}
\title{
Wine dataset
}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three different cultivars. 
The analysis determined the quantities of 13 constituents found in each of the three types of wines. 
}
\usage{data(wine)}
\format{
  A data frame with 178 observations on the following 14 variables :
  \describe{
    \item{\code{class}}{The class vector, the three different cultivars of wine are reprensented by the three integers : 1 to 3.}
	\item{\code{V1}}{Alcohol}
	\item{\code{V2}}{Malic acid}
	\item{\code{V3}}{Ash}
	\item{\code{V4}}{Alcalinity of ash}
	\item{\code{V5}}{Magnesium}
	\item{\code{V6}}{Total phenols}
	\item{\code{V7}}{Flavanoids}
	\item{\code{V8}}{Nonflavanoid phenols}
	\item{\code{V9}}{Proanthocyanins}
	\item{\code{V10}}{Color intensity}
	\item{\code{V11}}{Hue}
	\item{\code{V12}}{OD280/OD315 of diluted wines}
	\item{\code{V13}}{Proline}
  }
}

\source{
This dataset is from the UCI machine learning repository, provided here : \url{http://archive.ics.uci.edu/ml/datasets/Wine}.
}
\examples{
data(wine)
}
\keyword{datasets}
