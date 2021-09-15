\name{Crabs}
\alias{Crabs}
\docType{data}
\title{
Morphological Measurements on Leptograpsus Crabs.
}
\description{
The Crabs data frame has 200 rows and 6 columns, describing 5 morphological measurements on 50 crabs each of two colour forms and both sexes, of the species Leptograspus Variegatus collected at Fremantle, W. Australia.
}
\usage{data(Crabs)}
\format{
  A data frame with 200 observations on the following 6 variables.
  \describe{
    \item{\code{class}}{Type of the crabs: the first character represents the species - "B" or "O" for blue or orange-, the second represents the sex -"M" or "F" for male or female-.}
    \item{\code{FL}}{Frontal lob size (mm).}
    \item{\code{RW}}{Rear width (mm).}
    \item{\code{CL}}{Carapace length (mm).}
    \item{\code{CW}}{Carapace width (mm).}
    \item{\code{BD}}{Body depth (mm).}
  }
}
\details{
This dataset can also be found in the MASS package, the unique difference is the class vector which is easier to use here.
}
\source{
Campbell, N. A. and Mahon, R. J. (1974) \dQuote{A multivariate study of variation on two species of rock crab of genus Leptograspus}, \emph{Australian Journal of Zoology}, \bold{22}, 417--425.
}
\references{
Venables, W. N. and Ripley, B. D. (2002) \dQuote{Modern Applied Statistics with S}. Fourth edition. Springer.
}
\keyword{datasets}
