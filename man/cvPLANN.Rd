\name{cvPLANN}
\alias{cvPLANN}

\title{
Tune a Survival Neural Network Based on the PLANN Method
}

\description{
This function finds the optimal inter, size, decay, maxit, and MaxNWts parameters for the survival neural network by using cross-validation and the concordance index.
}

\usage{
cvPLANN(times, failures, pro.time=NULL, group=NULL, cov.quanti=NULL,
cov.quali=NULL, data, cv=10, inter, size, decay, maxit, MaxNWts)
}

\arguments{
  \item{times}{The name of the variable related the numeric vector with the follow-up times.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event).}
    \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "loglik", "ibs", "bll", and "ibll". Default value is the time at which half of the subjects are still at risk.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones. The default value is \code{NULL}: no specific exposure/treatment is considered. When a specific exposure/treatment is considered, it will be forced in the algorithm or related interactions will be tested when possible.}
  \item{cov.quanti}{The name(s) of the variable(s) related to the possible quantitative covariates. These variables must be numeric.}
  \item{cov.quali}{The name(s) of the variable(s) related to the possible qualitative covariates. These variables must be numeric with two levels: 0 and 1. A complete disjunctive form must be used for covariates with more levels.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time (\code{times}), the event (\code{failures}), the optional treatment/exposure (\code{group}) and the covariables included in the previous model (\code{cov.quanti} and \code{cov.quali}).}
  \item{cv}{The value of the number of folds. The default value is 10.}
  \item{inter}{The length of the intervals.}
  \item{size}{The number of units in the hidden layer.}
  \item{decay}{The parameter for weight decay.}
  \item{maxit}{The maximum number of iterations.}
  \item{MaxNWts}{The maximum allowable number of weights.}
 }


\value{
\item{optimal}{The value of \code{inter}, \code{size}, \code{decay}, \code{maxit}, and \code{MaxNWts} that gives the maximum mean cross-validated C-index.}
\item{results}{The data frame with the mean cross-validated C-index according to \code{inter}, \code{size}, \code{decay}, \code{maxit}, and \code{MaxNWts}.}
}

\references{
Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.
}

\examples{
data(dataK)

# The hyper-parameter grid needs to be more precise and the maximum number
# of iterations > 1000. We have reduced the arguments to respect examples requiring
# less than 5 seconds for packages on the CRAN.

tune.model <- cvPLANN(times="time", failures="event", data=dataK,
  cov.quali=c("sex", "stade", "delay"), cv=3,
  inter=365, size=c(16, 32), decay=0.01, maxit=50, MaxNWts=10000)

tune.model$optimal # the optimal hyperparameters

tune.model$results # the C-index for the tested grid
}

\keyword{Tunning}
