\name{rsPLANN}
\alias{rsPLANN}

\title{
Relative Survival Neural Network by Using the PLANN Method
}
\description{
This function allows to compute a predicting tool for relative survival additive risks with neural network using the PLANN method.
}

\usage{
rsPLANN(formula, data, pro.time=NULL, inter, size = 32, decay = 0.01,
                    maxit = 100, MaxNWts = 10000, trace = FALSE,
                    ratetable, age, year, sex)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function. The time must be in days.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{pro.time}{A numeric value with the maximum prognostic time in days. If \code{NULL}, the maximum prognostic time is the last time for wich an event is observed.}
\item{inter}{A numeric value representing the length of the intervals.}
\item{size}{A numeric value for the number of units in the hidden layer. Default is set to 32}
\item{decay}{A numeric value for the parameter for weight decay. Default is set to 0.01}
\item{maxit}{A numeric value for the maximum number of iterations. Default is set to 100.}
\item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the code, but increasing MaxNWts will probably allow fits that are very slow and time-consuming. Default is set to 10000}
\item{trace}{A logical argument for tracing optimization. Default is set to FALSE.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{A numeric vector with the patient age at baseline in days.}
\item{year}{A numeric vector with the baseline date in the study respecting the date format: number of days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
}

\value{  
\item{formula}{The \code{formula} object.}
\item{data}{The data frame given with the computed \code{Intervals} column added.}
\item{pro.time}{The maximum prognostic time.}
\item{fitsurvivalnet}{The fitted \code{nnet} object.}
\item{predictions}{A list with the predicted outcomes for the \code{N} individauls in \code{data}: \code{times} is a vector with the \code{P} prognostic times up to \code{pro.time}, \code{CIFc} is a \code{N x P} matrix with the cummulative probabilities of death caused by the disease, \code{FCIFp} is a \code{N x P} matrix with the cummulative probabilities of death of the reference population in \code{ratetable}, \code{tPcure} is a \code{N x P} matrix with the probabilities of being cured conditionnaly to be alive at the progostic times, \code{aPcure} is a vector with (asymptotic) probabilty of being cured, \code{Sc} is the survival among uncured patients (i.e. given the mortality is all releted to the disease), and Sp is the survival among cured patients (as in the reference population).}
}

\details{
The Biganzoli's PLANN method is used to predict the overall survival.
}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}


\examples{
data(dataK)
data(fr.ratetable)

cmp_fit <- cmp.rel(Surv(time, event) ~ 1, data=dataK[1:500,], ratetable=fr.ratetable,
        rmap=list(age=age, sex= sex, year=year))

plot(cmp_fit, col=c("black", "red"), lty=1, xscale = 365.241, xlab="Time (years)")

rs1 <- rsPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker, data = dataK[1:500,],
               pro.time = 365.241*20, inter=365.241, size = 16, decay = 0.01,
               maxit = 100, MaxNWts=10000, ratetable=fr.ratetable,
               age="age", sex="sexchara", year="year")

lines(rs1$times/365.241, rs1$mpredictions$CIF_C, type="s", col="black", lty=2)

lines(rs1$times/365.241, rs1$mpredictions$CIF_P, type="s", col="red", lty=2)

legend("topleft", legend=c("Cancer-related death", "Death due to other causes"), col=c("red", "blue"), lty=c(1,1))
}

\keyword{Neural Network}
