\name{cvPLANN}
\alias{cvPLANN}

\title{
Cross-validation method for a Neural Netwotk Model with the PLANN Method. 
}
\description{
Cross-validation method for the hyper-parameters for an object of class \code{sPLANN}.
}

\usage{
cvPLANN(formula, pro.time=NULL, data, cv=10, inter=1, size=8, decay=0.01,
                    maxit=1000, MaxNWts=10000)
}

\arguments{
  \item{formula}{The \code{formula} object.}
  \item{pro.time}{The prognostic time at which the metric is evaluated. If \code{NULL}, the median 
  of the times is chosen.}
  \item{data}{A data frame in which to look for the variables included in the \code{formula}.}
  \item{cv}{The number of splits for cross-validation. The default value is 10.}
  \item{inter}{A numeric value representing the length of the intervals.}
  \item{size}{A numeric value for the number of units in the hidden layer. Default is set to 8}
  \item{decay}{A numeric value for the parameter for weight decay. Default is set to 0.01}
  \item{maxit}{A numeric value for the maximum number of iterations. Default is set to 1000.}
  \item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the code, but increasing MaxNWts will probably allow fits that are very slow and time-consuming. Default is set to 10000}
}

\value{
  \item{optimal}{A list giving the optimal value for each parameter.}
  \item{results}{A data frame listing the parameters combinaison and their metrics values.}
}

\examples{
data(dataK) # the database with the observed sample

tune.sPLANN <- cvPLANN(Surv(time, event)~ stade + delay, data=dataK, cv=3,
                     inter=365.241, size=c(2, 4), decay=0.01)
                     
tune.sPLANN$optimal
}

\seealso{ \code{\link{sPLANN}} }

\keyword{Cross-validation}
