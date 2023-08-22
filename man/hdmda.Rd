\name{hdmda}
\alias{hdmda}
\title{Mixture Discriminant Analysis with HD Gaussians}
\description{HD-MDA implements mixture discriminant analysis (MDA, Hastie & Tibshirani, 1996) with HD Gaussians instead of full Gaussians. Each class is assumed to be made of several class-specific groups in which the data live in low-dimensional subspaces. From a technical point of view, a clustering is done using \code{\link{hddc}} in each class.}
\usage{
hdmda(X,cls,K=1:10,model='AkjBkQkDk',show=FALSE,...)
}
\arguments{
  \item{X}{A matrix or a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.}
   \item{cls}{The vector of the class of each observations, its type can be numeric or string.}
   \item{K}{A vector of integers specifying the number of clusters for which the BIC and the parameters are to be calculated; the function keeps the parameters which maximises the BIC. Note that the length of the vector K can't be larger than 20. Default is 1:10.}
   \item{model}{A character string vector, or an integer vector indicating the models to be used. The available models are: "AkjBkQkDk" (default), "AkBkQkDk", "ABkQkDk", "AkjBQkDk", "AkBQkDk", "ABQkDk", "AkjBkQkD", "AkBkQkD", "ABkQkD", "AkjBQkD", "AkBQkD", "ABQkD", "AjBQD", "ABQD". It is not case sensitive and integers can be used instead of names, see details for more information. Several models can be used, if it is, only the results of the one which maximizes the BIC criterion is kept. To run all models, use model="ALL".}
   \item{show}{Use show = TRUE to display some information related to the clustering.}
   \item{...}{Any argument that can be used by the function \code{\link{hddc}}.}
}
\value{
\code{hdmda} returns an 'hdmda' object which is a list containing:
\item{alpha}{Estimated prior probabilities for the classes.}
\item{prms}{Estimated mixture parameters for each class.}
\item{kname}{The name (level) of each class.}
}
\references{

C. Bouveyron and C. Brunet (2014), \dQuote{Model-based clustering of high-dimensional data: A review}, Computational Statistics and Data Analysis, vol. 71, pp. 52-78.

Bouveyron, C. Girard, S. and Schmid, C. (2007), \dQuote{High Dimensional Discriminant Analysis}, Communications in Statistics: Theory and Methods, vol. 36 (14), pp. 2607-2623.

Bouveyron, C. Celeux, G. and Girard, S. (2011), \dQuote{Intrinsic dimension estimation by maximum likelihood in probabilistic PCA}, Pattern Recognition Letters, vol. 32 (14), pp. 1706-1713.

Berge, L. Bouveyron, C. and Girard, S. (2012), \dQuote{HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data}, Journal of Statistical Software, 46(6), pp. 1-29, url: \doi{10.18637/jss.v046.i06}.

Hastie, T., & Tibshirani, R. (1996), \dQuote{Discriminant analysis by Gaussian mixtures}, Journal of the Royal Statistical Society, Series B (Methodological), pp. 155-176.
}

\details{
Some information on the signification of the model names:
\describe{
	\item{Akj are the parameters of the classes subspaces:}{
		\itemize{
		\item{if Akj: each class has its parameters and there is one parameter for each dimension}
		\item{if Ak: the classes have different parameters but there is only one per class}
		\item{if Aj: all the classes have the same parameters for each dimension (it's a particular case with a common orientation matrix)}
		\item{if A: all classes have the same one parameter}
		}
	}

	\item{Bk are the noises of the classes subspaces:}{
		\itemize{
			\item{If Bk: each class has its proper noise}
			\item{if B:  all classes have the same noise}
		}
	}

	\item{Qk is the orientation matrix of each class:}{ 
		\itemize{
			\item{if Qk: all classes have its proper orientation matrix}
			\item{if Q: all classes have the same orientation matrix}
		}
	}
	
	\item{Dk is the intrinsic dimension of each class:}{ 
		\itemize{
			\item{if Dk: the dimensions are free and proper to each class}
			\item{if D: the dimension is common to all classes}
		}
	}
}
The model \dQuote{all} will compute all the models, give their BIC and keep the model with the highest BIC value.
Instead of writing the model names, they can also be specified using an integer.  1 represents the most general model (\dQuote{AkjBkQkDk}) while 14 is the most constrained (\dQuote{ABQD}), the others  number/name matching are given below:
\tabular{lcclc}{
AkjBkQkDk \tab   1   \tab   \tab  AkjBkQkD \tab   7   \cr 
AkBkQkDk \tab   2   \tab \tab  AkBkQkD \tab   8   \cr   
ABkQkDk \tab   3   \tab  \tab ABkQkD \tab   9   \cr   
AkjBQkDk \tab   4   \tab  \tab  AkjBQkD \tab   10   \cr   
AkBQkDk \tab   5   \tab  \tab  AkBQkD \tab   11   \cr   
ABQkDk \tab   6   \tab  \tab  ABQkD \tab   12  \cr
AjBQD \tab 13 \tab  \tab ABQD \tab 14
}

}

\author{
Laurent Berge, Charles Bouveyron and Stephane Girard
}

\seealso{\code{\link{hdda}}, \code{\link{hddc}}}
\examples{
# Load the Wine data set
data(wine)
cls = wine[,1]; X = scale(wine[,-1])

# A simple use...
out = hdmda(X[1:100,],cls[1:100])
res = predict(out,X[101:nrow(X),])

# Comparison between hdmda and hdda in a CV setup
set.seed(123); nb = 10; Err = matrix(NA,2,nb)
for (i in 1:nb){
  cat('.')
  test = sample(nrow(X),50)
  out0 = lda(X[-test,],cls[-test])
  res0 = predict(out0,X[test,])
  Err[1,i] = sum(res0$class != cls[test]) / length(test)
  out = hdmda(X[-test,],cls[-test],K=1:3,model="AKJBQKDK")
  res = predict(out,X[test,])
  Err[2,i] = sum(res$class != cls[test]) / length(test)
}
cat('\n')
boxplot(t(Err),names=c('LDA','HD-MDA'),col=2:3,ylab="CV classifciation error",
  main='CV classifciation error on Wine data')
}

\keyword{mixture discriminant analysis}
\keyword{high-dimensional data}
