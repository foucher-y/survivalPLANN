
cvPLANN <- function(times, failures, pro.time=NULL, group=NULL, cov.quanti=NULL, cov.quali=NULL,
                           data, cv=10, inter, size, decay, maxit, MaxNWts){ # @Thomas -> Passer cette fonction avec Surv(...) 

  data.plann <- data[,c(times, failures, group, cov.quanti, cov.quali)]
  
  if(is.null(pro.time)) {pro.time <- median(data[,times])}
  
  sample_id <- sample(nrow(data.plann))
  folds <- cut(seq(1,nrow(data.plann)), breaks=cv, labels=FALSE)
  folds_id <- folds[sample_id]
  data.plann$folds <- folds_id

  #.f  <- as.formula(paste("Surv(", times, ",", failures, ")", "~."))

  .outcome <- paste("Surv(", times, ",", failures, ")")

  if(!(is.null(group))) {
    if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", group, "+", paste( cov.quanti,  collapse = " + "),
                              " + ", paste(cov.quali, collapse = " + "),
                              collapse = " ") )
    }
    if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", group, "+",
                              paste( cov.quanti, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", group, "+",
                              paste(cov.quali, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", group) )
    }
  }   else {
    if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", paste( cov.quanti,  collapse = " + "), " + ", paste(cov.quali, collapse = " + "),
                              collapse = " ") )
    }
    if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", paste( cov.quanti, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~",  paste(cov.quali, collapse = " + "),collapse = " ") )
    }
  }


  .time <- sort(unique(data.plann[,times]))

  .grid <-  expand.grid(inter=inter, size=size, decay=decay, maxit=maxit, MaxNWts=MaxNWts)

  .CVtune<-vector("list",cv*dim(.grid)[1])

  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .CVtune[[l]]<-list(train=data.plann[data.plann$folds!=k, ], valid=data.plann[data.plann$folds==k, ], grid=.grid[j,])
      l=l+1
    }
  }

  plann.time.par<-function(xx, times, failures, group, cov.quanti, cov.quali, newtimes){

    inter=xx$grid$inter
    size=xx$grid$size
    decay=xx$grid$decay
    maxit=xx$grid$maxit
    MaxNWts=xx$grid$MaxNWts

    data=xx$train
    newdata=xx$valid

    if(!(is.null(group))){
      .data <- data[,c(times, failures, group, cov.quanti, cov.quali)]}   else{
        .data <- data[,c(times, failures, cov.quanti, cov.quali)] }

    #.f  <- as.formula(paste("Surv(", times, ",", failures, ")", "~."))

    .outcome <- paste("Surv(", times, ",", failures, ")")

    if(!(is.null(group))) {
      if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
        .f <- as.formula( paste(.outcome, "~", group, "+", paste( cov.quanti,  collapse = " + "),
                                " + ", paste(cov.quali, collapse = " + "),
                                collapse = " ") )
      }
      if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
        .f <- as.formula( paste(.outcome, "~", group, "+",
                                paste( cov.quanti, collapse = " + "),collapse = " ") )
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
        .f <- as.formula( paste(.outcome, "~", group, "+",
                                paste(cov.quali, collapse = " + "),collapse = " ") )
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==T) {
        .f <- as.formula( paste(.outcome, "~", group) )
      }
    }   else {
      if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
        .f <- as.formula( paste(.outcome, "~", paste( cov.quanti,  collapse = " + "), " + ", paste(cov.quali, collapse = " + "),
                                collapse = " ") )
      }
      if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
        .f <- as.formula( paste(.outcome, "~", paste( cov.quanti, collapse = " + "),collapse = " ") )
      }
      if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
        .f <- as.formula( paste(.outcome, "~",  paste(cov.quali, collapse = " + "),collapse = " ") )
      }
    }

    .plann <- sPLANN(.f, data = .data,
                  inter = inter, size = size, decay = decay,  maxit = maxit, MaxNWts = MaxNWts)

    .time<-sort(unique(.data[,times]))

    .newdata <- data.frame(newdata[,c(group, cov.quanti, cov.quali)])
    .pred.temp <- predict(.plann, newdata=newdata)
    .pred <- .pred.temp$predictions
    .time.plann <- .pred.temp$times

    if(!is.null(newtimes)) {
      .pred.plann <- cbind(rep(1, dim(.pred)[1]), .pred)
      .time.plann <- c(0, .time.plann)

      idx=findInterval(newtimes, .time.plann)
      .pred=.pred.plann[,pmax(1,idx)]

      .time <- newtimes
    }

    return(as.matrix(.pred))
  }

  .preFIT<-list()
  .preFIT<-lapply(.CVtune, plann.time.par, times=times, failures=failures, group=group,
                 cov.quanti=cov.quanti, cov.quali=cov.quali, newtimes=.time)

  .FitCV <- replicate(dim(.grid)[1], matrix(NA, nrow = length(data[,times]),
                                            ncol = length(.time)), simplify=F)
  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .FitCV[[j]][data.plann$folds==k,] <- .preFIT[[l]]
      l<-l+1
    }
  }

  plann.best.measure <- function(prediction.matrix, times, failures, data, prediction.times){
    return(metrics(times=times, failures=failures, data=data, prediction.matrix=prediction.matrix,
                  prediction.times=prediction.times, pro.time=pro.time, metric="ci"))
  }

    .measure<-sapply(.FitCV, plann.best.measure, times=times, failures=failures, data=data.plann, prediction.times=.time)

    .res <- data.frame(inter = .grid[,1], size = .grid[,2], decay=.grid[,3],
                       maxit = .grid[,4], MaxNWts = .grid[,5], measure = .measure)

    .maxi<-.res[which(.res$measure==max(.res$measure, na.rm=TRUE) & is.na(.res$measure)==FALSE),]
    .maxi<-.maxi[1,]

    return( list(optimal=list(inter=.maxi$inter,
                              size=.maxi$size,
                              decay=.maxi$decay,
                              maxit=.maxi$maxit,
                              MaxNWts=.maxi$MaxNWts),
                              results=.res ))
}
