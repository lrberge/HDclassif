\name{HDclassif-package}
\alias{HDclassif-package}
\alias{HDclassif}
\docType{package}
\title{
High Dimensional Discriminant Analysis and Data Clustering
}
\description{
Disciminant analysis and data clustering methods for high dimensional data, based on the asumption that high-dimensional data live in different subspaces with low dimensionality, proposing a new parametrization of the Gaussian mixture model which combines the ideas of dimension reduction and constraints on the model.
}
\details{
\tabular{ll}{
Package: \tab HDclassif\cr
Type: \tab Package\cr
Version: \tab 2.1.0\cr
Date: \tab 2018-05-11\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package is used to make efficient supervised and unsupervised classification with high dimensional data.
The supervised method uses the \var{hdda} function to get the data parameters and the \var{predict} function to realise the class prediction of a dataset.
The unsupervised method is implemented in the \var{hddc} function, and once the parameters are estimated, the \var{predict} gives the class prediction of other datasets. The method used in the hddc is based on the Expectation - Maximisation algorithm.
}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard

Maintainer: Laurent Berge <laurent.berge at uni.lu> 
}
\references{
Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High Dimensional Discriminant Analysis}, \emph{Communications in Statistics: Theory and Methods}, vol. \bold{36} (14), pp. 2607--2623

Bouveyron, C. Girard, S. and Schmid, C. (2007) \dQuote{High-Dimensional Data Clustering}, \emph{Computational Statistics and Data Analysis}, vol. \bold{52} (1), pp. 502--519

Berge, L. Bouveyron, C. and Girard, S. (2012) \dQuote{HDclassif: An R Package for Model-Based 
Clustering and Discriminant Analysis of High-Dimensional Data}, \emph{Journal of Statistical Software}, 
\bold{46}(6), 1--29, url: \href{https://doi.org/10.18637/jss.v046.i06}{https://doi.org/10.18637/jss.v046.i06}
}

\keyword{ package }
