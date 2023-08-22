\name{predict.hdmda}
\alias{predict.hdmda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction method for \sQuote{hdmda} class objects.
}
\description{
This function computes the class prediction of a dataset with respect to the model-based supervised classification method \code{\link{hdmda}}.
}
\usage{
\method{predict}{hdmda}(object, X, \dots)
}

\arguments{
  \item{object}{
An object of class \sQuote{hdmda}.
}
  \item{X}{
A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{\dots}{
Arguments based from or to other methods. Not currently used.
}
}

\value{
	\item{class}{vector of the predicted class.}
	\item{posterior}{The matrix of the probabilities to belong to a class for each observation and each class.}
}
\references{
	C. Bouveyron and C. Brunet (2014), \dQuote{Model-based clustering of high-dimensional data: A review}, Computational Statistics and Data Analysis, vol. 71, pp. 52-78.
	
	Bouveyron, C. Girard, S. and Schmid, C. (2007), \dQuote{High Dimensional Discriminant Analysis}, Communications in Statistics: Theory and Methods, vol. 36 (14), pp. 2607-2623.
	
	Bouveyron, C. Celeux, G. and Girard, S. (2011), \dQuote{Intrinsic dimension estimation by maximum likelihood in probabilistic PCA}, Pattern Recognition Letters, vol. 32 (14), pp. 1706-1713.
	
	Berge, L. Bouveyron, C. and Girard, S. (2012), \dQuote{HDclassif: An R Package 
	for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, 
	Journal of Statistical Software, 46(6), pp. 1-29, url: \doi{10.18637/jss.v046.i06}.
	
	Hastie, T., & Tibshirani, R. (1996), \dQuote{Discriminant analysis by Gaussian mixtures}, Journal of the Royal Statistical Society, Series B (Methodological), pp. 155-176.
}

\author{
Laurent Berge, Charles Bouveyron and Stephane Girard
}

\seealso{
\code{\link{hdmda}}
}
\examples{
# Load the Wine data set
data(wine)
cls = wine[,1]; X = scale(wine[,-1])

# A simple use...
out = hdmda(X[1:100,],cls[1:100])
res = predict(out,X[101:nrow(X),])
}
