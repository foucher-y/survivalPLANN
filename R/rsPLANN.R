
rsPLANN <- function(formula, data, pro.time=NULL, inter, size= 32, decay=0.01,
                    maxit=100, MaxNWts=10000, trace=FALSE,
                    ratetable, age, year, sex) # @Thomas: inclure ratetable et ses arguments dans formula
{
  
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(inter)) stop("an inter argument is required")
  if (as.character(class(formula)) != "formula") stop("The formula argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (as.character(class(inter)) != "numeric") stop("The inter argument must be numeric")
  
  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
  ####### data management
  
  splann <- sPLANN(formula, data=data, pro.time=pro.time, inter=inter, 
                          size=size, decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  predO <- predict(splann, newtimes=splann$intervals)
  
  times <- predO$times
  
  #expectedhaz <- function(ratetable, age, year, sex, time) 
  #{
  #  time <- min(time, 1000000)
  #  .year <- date.mdy(year+time)$year
  #  
  #  ratetable[as.character( min( floor((age+time)/365.24), max(as.numeric(names(ratetable[, "2000", "male"]))) ) ),
  #            as.character( min( .year, max(as.numeric(names(ratetable["51", , "male"]))) ) ),
  #            sex]
  #}
  
  exphaz <- function(x,age,sex,year) { survivalNET::expectedhaz(ratetable, age=age, sex=sex, year=year, time=x)}
  
  survO <- as.matrix(predO$predictions)
  dimnames(survO) <- NULL
  
  N <- dim(survO)[1]
  P <- dim(survO)[2]
  
  hP <- matrix(-99, ncol=length(times), nrow=N)
  
  for (i in 1:N) # @Thomas : merci de voir si tu augmenter la vitesse du calcul de hP
  {
    hP[i,] <- sapply(times, FUN="exphaz", age=data[i,age],
                     sex=data[i,sex], year=data[i,year]) * splann$inter
  }
  
  hcumO <- -1*log(survO)
  hinstO <- hcumO[,2:length(times)] - hcumO[,1:(length(times)-1)]
  hinstO[hinstO==Inf] <- NA
  
  for (i in 1:N)
  {
    if(sum(survO[i,]==0)>0)
    {
      hinstO[i,is.na(hinstO[i,])] <- hinstO[i,!is.na(hinstO[i,])][sum(!is.na(hinstO[i,]))]
    }
  }
  
  distOa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstO )))))
  distOb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstO )))))
  distO <- cbind(rep(0, N), (distOa + distOb)/2)
  
  hinstP <- hP[,1:(length(times)-1)]
  distPa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstP )))))
  distPb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstP )))))
  distP <- cbind(rep(0, N), (distPa + distPb)/2)
  
  hinstE <- hinstO - hinstP
  distEa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstE )))))
  distEb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstE )))))
  distE <- cbind(rep(0, N), (distEa + distEb)/2)
  
  distP <- distP * (1-survO)/distO
  distE <- distE * (1-survO)/distO
  
  distP[survO==1] <- 0
  distE[survO==1] <- 0
  
  distO <- distP + distE
  
  distPinf <- distP[,P]
  distEinf <- distE[,P]
  
  estimPcure <- (round(distPinf + distEinf, 10) == 1)
  
  survP <- cbind(rep(1, N), exp(-t(as.matrix(cumsum(data.frame(t(hinstP)))))))
  survU <- cbind(rep(1, N), exp(-t(as.matrix(cumsum(data.frame(t(hinstE)))))))
  
  Pcure <- distPinf / (distPinf + (1-distPinf) * survU)

  sumS <- apply((1-distO), FUN="sum", MARGIN=2)
  
  SlE <- (1-distO) * cbind(rep(0, N), hinstE)
  sumSlE <- apply(SlE, FUN="sum", MARGIN=2)
  lambda_C <- sumSlE/sumS
  #Lambda_C <- cumsum(lambda_C)
  
  SlP <- (1-distO) * cbind(rep(0, N), hinstP)
  sumSlP <- apply(SlP, FUN="sum", MARGIN=2)
  lambda_P <- sumSlP/sumS
  #Lambda_P <- cumsum(lambda_P)
  
  # warning -> NA pour tCure ...
  
  res <- list(formula = formula,
              data = data,
              ratetable = ratetable,
              age = age,
              sex= sex,
              year = year,
              pro.time = pro.time,
              inter = splann$inter,
              size = splann$size,
              decay = splann$decay,
              fitsurvivalnet = splann,
              times = times,
              ipredictions = list(survival_P=survP,
                                  survival_O=1-distO,
                                  survival_E2=(1-distO)/survP,
                                  survival_E=survU, # remarque : S(1-distO)/survP = survU
                                  CIF_C = distE, CIF_P = distP, maxCIF_P = distPinf,
                                  cure = Pcure),
              mpredictions = list(survival_O = apply((1-distO), FUN="mean", MARGIN=2),
                                  survival_P = apply(survP, FUN="mean", MARGIN=2),
                                  survival_R = apply((1-distO), FUN="mean", MARGIN=2)/apply(survP, FUN="mean", MARGIN=2),
                                  survival_E = apply((1-distO)/survP, FUN="mean", MARGIN=2),
                                  CIF_C =  apply(distE, FUN="mean", MARGIN=2),
                                  CIF_P =  apply(distP, FUN="mean", MARGIN=2),   
                                  cure = apply(Pcure, FUN="mean", MARGIN=2),
                                  hazard_P = lambda_P,
                                  hazard_C = lambda_C
                                  )
              )
  
  class(res) <- "rsPLANN"
  return(res)
}
