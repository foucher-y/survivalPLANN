\name{predict.rsPLANN}
\alias{predict.rsPLANN}

\title{
Predictions From a Relative Neural Network Based on the PLANN Method
}

\description{
This function produces several predictions from a relative neural network based on the PLANN method. 
}

\usage{
\method{predict}{rsPLANN}(object, newdata = NULL, newtimes = NULL, ratetable, age, year, sex, ...)
}

\arguments{
\item{object}{The result of the \code{rsPLANN} function.}
\item{newdata}{An optional data frame comprising of new examples to be predicted. If \code{NULL}, the data frame used is the one used for training in the \code{rsPLANN} function.}
\item{newtimes}{A optional numeric vector comprising of times to get survival estimations. If \code{NULL}, the times are the intervals used in the \code{rsPLANN} function.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{A numeric vector with the patient age at baseline in days.}
\item{year}{A numeric vector with the baseline date in the study respecting the date format: number of days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
\item{...}{Further arguments passed to or from other methods.}
}


\value{
\item{times}{The \code{P} times used for the predicitions.}
\item{predictions}{A list with the predicted outcomes for the \code{N} individauls: \code{Fc} is a \code{N x P} matrix with the cummulative probabilities of death caused by the disease, \code{Fp} is a \code{N x P} matrix with the cummulative probabilities of death of the reference population in \code{ratetable}, \code{tPcure} is a \code{N x P} matrix with the probabilities of being cured conditionnaly to be alive at the progostic times, and \code{aPcure} is a vector with (asymptotic) probabilty of being cured.}
}


\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}

\examples{
data(dataK)
data(fr.ratetable)

rs1 <- rsPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker, data = dataK[1:300,],
               pro.time = 365.241*20, inter=365.241/4, size = 16, decay = 0.01,
               maxit = 100, MaxNWts=10000, ratetable=fr.ratetable,
               age="age", sex="sexchara", year="year")
               
p <- predict(rs1, newdata=dataK[301:500,],
   ratetable=fr.ratetable,  age="age", sex="sexchara", year="year")
   
plot(rs1$times/365.241, rs1$mpredictions$CIF_C,
     xlab="Post-diagnostic time (years)",
     ylab="Cumulative incidence functions",
     ylim=c(0,1), col="black", type="s",)

lines(rs1$times/365.241, rs1$mpredictions$CIF_P, type="s", col="red", lty=2)

legend("topleft", legend=c("Cancer-related death", "Death due to other causes"),
       col=c("red", "blue"), lty=c(1,1))
}

\keyword{Predicting}
