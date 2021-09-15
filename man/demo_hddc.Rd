\name{demo_hddc}
\alias{demo_hddc}
\title{
Demonstration of the clustering process of HDDC.
}
\description{
This demonstration uses a PCA on the first two principal axis of the Crabs dataset -that can be found in the package- to show the clustering process of HDDC. At each step of the clustering, the means and directions are shown by, respectively, points and lines.
This function should only be used in demo(hddc).
}
\usage{
demo_hddc()
}
\note{
The algorithm and the initialization are interactively chosen.
}
\value{ The plots of the clustering process.}
\author{
Laurent Berge, Charles Bouveyron and Stephane Girard 
}
\seealso{
\code{\link{hddc}}.
}
\keyword{ demo }
