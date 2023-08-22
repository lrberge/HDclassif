\name{plot.hdc}
\alias{plot.hdc}
\title{
Cattell's Scree-Test for 'hdc' class objects.
}
\description{
This function plots Cattell's scree-test or the BIC selection, using parameters coming from \var{hdda} or \var{hddc} functions.
}
\usage{
\method{plot}{hdc}(x, method = NULL, threshold = NULL, noise.ctrl=1e-8, \dots)
}

\arguments{
  \item{x}{
A 'hdc' class object obtained using \var{hdda} or \var{hddc} methods.
}
	\item{method}{
	The method used to select the intrinsic dimension. It can be "BIC" or "Cattell". By default it takes the method used when obtaining the parameters using \var{hdda} or \var{hddc}.
	Note that "Cattell" (resp. "BIC") can be abreviated to "C" (resp. "B") and that this argument is not case sensitive.
}
  \item{threshold}{
The threshold used in Cattell's Scree-Test. By default it takes the threshold in the argument x, if none, the default value of the threshold is 0.2.
}
	\item{noise.ctrl}{
	This parameter avoids to have a too low value of the 'noise' parameter b. It garantees that the dimension selection process do not select too many dimensions (which leads to a potential too low value of the noise parameter b). When selecting the intrinsic dimensions using Cattell's scree-test or BIC, the function doesn't use the eigenvalues inferior to noise.ctrl, so that the intrinsic dimensions selected can't be higher or equal to the order of these eigenvalues.
}
\item{\dots}{
Arguments based from or to other methods.
}
}
\value{
\describe{
	\item{If method = "Cattell"}{ The plot of the eigen values and of the sequential differences of the eigen values. The dimension to retain is the one before the last fall of the eigenvalues' differences below the threshold.}
	\item{If method = "BIC"}{The BIC related to the dimension for each class. It stops after the first fall of the BIC.}
}
}
\references{
Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607--2623

Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519

Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based 
Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, 
\bold{46}(6), 1--29, url: \doi{10.18637/jss.v046.i06}
}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard
}

\seealso{
\code{\link{hdda}}, \code{\link{hddc}}, \code{\link{predict.hdc}}.
}
\examples{
# Example 1 :
data(wine)
a <- wine[,-1]
z <- wine[,1]

prms1 <- hdda(a, z, model="AkBkQkDk", scaling=TRUE, d_select="bic")

#the plot related to the selection that has been done: BIC
plot(prms1)     

#it shows the plot of Cattell's scree-test, with a threshold of .3
plot(prms1,"Cattell",0.3)                         


prms2 <- hdda(a, z, model="AkBkQkD", scaling=TRUE, d_select="cattell")
#the plot related to the selection that has been done: Cattell's scree-test
plot(prms2) 
#the plot of the BIC
plot(prms2,"b") 

}
\keyword{ cattell }
\keyword{ hdda }
\keyword{ hddc }
\keyword{ clustering }
