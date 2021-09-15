\name{simuldata}
\alias{simuldata}
\title{
Gaussian Data Generation
}
\description{
This function generates two datasets according to the model [AkBkQkDk] of the HDDA gaussian mixture model paramatrisation (see ref.).
}
\usage{
simuldata(nlearn, ntest, p, K = 3, prop = NULL, d = NULL, a = NULL, b = NULL)
}
\arguments{
  \item{nlearn}{
The size of the learning dataset to be generated.
}
  \item{ntest}{
The size of the testing dataset to be generated.
}
  \item{p}{
The number of variables.
}
  \item{K}{
The number of classes.
}
  \item{prop}{
The proportion of each class.
}
  \item{d}{
The dimension of the intrinsic subspace of each class.
}
  \item{a}{
The value of the main parameter of each class.
}
  \item{b}{
The noise of each class.
}
}
\value{
\item{ X }{ The learning dataset.}
\item{ clx }{ The class vector of the learning dataset. }
\item{ Y }{ The test dataset. }
\item{ cly }{ The class vector of the test dataset. }
\item{ prms }{ The principal parameters used to generate the datasets. }
}
\references{
Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics : Theory and Methods}, vol. \bold{36}(14), pp. 2607--2623
}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard 
}
\seealso{
\code{\link{hddc}}, \code{\link{hdda}}.
}
\examples{
data <- simuldata(500, 1000, 50, K=5, prop=c(0.2,0.25,0.25,0.15,0.15))
X <- data$X
clx <- data$clx
f <- hdda(X, clx)
Y <- data$Y
cly <- data$cly
e <- predict(f, Y, cly)

}
\keyword{ generation }
\keyword{ gaussian }
