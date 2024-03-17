
require(plyr)

bestbf <- function(BF,digits=2,restrict=NA,save=FALSE) {
  if (!save) cat("BEST MODEL\n")
  outs <- list(best=NULL,models=NULL)
  if ( max(as.vector(BF))<1 ) {
    out <- sort(1/as.vector(BF),decreasing=F)
    if (save) {
      outs$best=c(Intercept=1)
      if (is.na(restrict)) outs$models <- round(out,digits=digits) else
        outs$models <- round(out[out<restrict],digits=digits)

    } else {
      cat("Intercept\n")
      if (is.na(restrict)) print(round(out,digits=digits)) else
        print(round(out[out<restrict],digits=digits))
    }
  } else {
    out <- sort(as.vector(BF[which.max(BF)])/as.vector(BF),
                     decreasing=F)
    if (save) {
      outs$best <- round(as.vector(BF[which.max(BF)]),digits=digits) 
      if (is.na(restrict)) outs$models <- round(out,digits=digits) else
        outs$models <- round(out[out<restrict],digits=digits)      
    } else {
      print(round(as.vector(BF[which.max(BF)]),digits=digits)) 
      out <- sort(as.vector(BF[which.max(BF)])/as.vector(BF),
                     decreasing=F)
      if (is.na(restrict)) print(round(out,digits=digits)) else
        print(round(out[out<restrict],digits=digits))
    }
  }
  invisible(outs)
}

bma <- function(bf) {

  logSumExpLogs <- function(v)
    logMeanExpLogs(v) + log(length(v))

  logPostProbs <- function(bf){
    lbfs = c(extractBF(bf,logbf=TRUE)$bf,0)
    lpps = lbfs - logSumExpLogs(lbfs)
    names(lpps) = c(names(bf)$numerator,"Intercept")
    return(lpps)
  }

  bmaf <- function(bf,f) {
    mnames <- names(bf)$numerator
    pps <- exp(logPostProbs(bf))
    isin <- c(c(1:length(mnames)) %in% grep(f,mnames),FALSE)
    (sum(pps[isin])/sum(pps[!isin]))/(sum(isin)/sum(!isin))
  }
  
  terms <- unique(unlist(strsplit(names(bf)[[1]]," + ",fixed=T)))
  terms <- terms[order(unlist(lapply(strsplit(terms,":",fixed=TRUE),length)))]
  out <- vector(mode="list",length=length(terms))
  names(out) <- terms
  for (i in terms) {
    out[[i]] <- bmaf(bf,i)
  }
  out <- unlist(out)
  BFfor <- out[out>=1]
  BFagainst <- 1/out[out<1]
  cat("BAYES FACTORS in favor of term INCLUSION\n")
  print(sort(BFfor,decreasing = TRUE))
  cat("\nBAYES FACTORS in favor of term EXCLUSION\n")
  print(sort(BFagainst,decreasing = FALSE))
  invisible(out)
}

library(car)

wsAnova=function(dat,SStype=3,spss=F) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  dat=dat[do.call(order,dat[,-dim(dat)[2]]),]
  snams=levels(dat[,1]); ns=length(snams)
  dvnam=names(dat)[dim(dat)[2]]
  facn=names(dat)[-c(1,dim(dat)[2])]
  nifac=length(facn)
  idata=data.frame(dat[dat[,1]==snams[1],facn])
  names(idata)=facn
  for (i in facn) 
    if (i==facn[1]) ifacn=as.character(idata[,1]) else 
      ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
  facnr=facn[nifac:1]
  e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
              ncol=length(ifacn),dimnames=list(snams,ifacn))
  print(summary(Anova(lm(e.mv ~ 1),
                idata=idata,type=SStype, 
                idesign=formula(paste("~",paste(facn,collapse="*")))),
          multivariate=FALSE))
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}

toSpss=function(dat) {
  if (dim(dat)[2]==2) return(dat)
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  dat=dat[do.call(order,dat[,-dim(dat)[2]]),]
  snams=levels(dat[,1]); ns=length(snams)
  dvnam=names(dat)[dim(dat)[2]]
  facn=names(dat)[-c(1,dim(dat)[2])]
  nifac=length(facn)
  idata=data.frame(dat[dat[,1]==snams[1],facn])
  names(idata)=facn
  for (i in facn) 
    if (i==facn[1]) ifacn=as.character(idata[,1]) else 
      ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
  facnr=facn[nifac:1]
  e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
              ncol=length(ifacn),dimnames=list(snams,ifacn))
  e.mv=cbind.data.frame(s=row.names(e.mv),e.mv)
  row.names(e.mv)=NULL
  e.mv
}



mneffects <- function(df,elist,digits=3,dvnam="y",save=FALSE,verbose=TRUE,
                   vars=F,err=F,log=F) {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam <- dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  for (i in 1:length(elist)) {
    if (verbose) cat(paste(paste(elist[[i]],collapse=":"),"\n"))
    mns <- tapply(df[,dvnam],df[,elist[[i]]],mean,na.rm=T)
    if (err) mns <- pnorm(mns)
    if (vars) mns <- sqrt(mns)
    if (log) mns <- exp(mns)
    if (verbose) {
      print(round(mns,digits))    
      cat("\n")
    }
  } 
  if (save) mns
}

sdeffects <- function(df,elist,digits=3,dvnam="y",save=FALSE,verbose=TRUE,
                      vars=F,err=F,log=F) {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam <- dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  for (i in 1:length(elist)) {
    if (verbose) cat(paste(paste(elist[[i]],collapse=":"),"\n"))
    sds <- tapply(df[,dvnam],df[,elist[[i]]],sd,na.rm=T)
    if (err) sds <- pnorm(sds)
    if (vars) sds <- sqrt(sds)
    if (log) sds <- exp(sds)
    if (verbose) {
      print(round(sds,digits))    
      cat("\n")
    }
  } 
  if (save) sds
}

se=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE",err=TRUE) {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  mn=tapply(df[,dvnam],df[,facs],mean)
  se=tapply(df[,dvnam],df[,facs],sd)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (is.na(ci)) mn else {
    if (ci=="SE") se/sqrt(ns) else
     qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
  }
}

add.bars=function(mn,se,xvals=NA,len=.1,antiprobit=FALSE,col="black") {
  
  plotbars <- function(x,m,l,h,len,col="black") {
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],l[j],length=len,angle=90,col=col)
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],h[j],length=len,angle=90,col=col)    
  }
  
  if (any(is.na(xvals))) if (is.matrix(mn)) 
    xvals <- as.numeric(dimnames(mn)[[2]]) else
    xvals <- as.numeric(names(mn))
  lo <- mn-se
  hi <- mn+se
  if (antiprobit) {
    mn=pnorm(mn)
    lo=pnorm(lo)
    hi=pnorm(hi)
  }
  if (!is.matrix(mn)) 
      plotbars(xvals,mn,lo,hi,col=col,len=len) else
    for (i in 1:dim(mn)[1]) 
      plotbars(x=xvals,m=mn[i,],l=lo[i,],h=hi[i,],len=len,col=col)
}    

arr2df=function(arr) {
  if (is.null(dim(arr))) out=data.frame(y=arr) else {
    dn=dimnames(arr)
    if (length(dn)==1) {
      out=cbind.data.frame(factor(dn[[1]],dn[[1]]),arr)
      names(out)=c(names(dn),"y")
      row.names(out)=NULL
    } else {
      tmp=vector(mode="list",length=length(dn))
      names(tmp)=names(dn)
      k=1
      for (j in names(dn)) {
        n=length(dn[[j]])
        tmp[[j]]=gl(n,k,length(arr),dn[[j]])
        k=k*n
      }
      out=cbind(data.frame(tmp),y=as.vector(arr))
      row.names(out)=NULL
    }
  }
  out
}

# spss=F; dvnam=NA; SStype=3; sfac="s"
mixedAnova=function(dat,bsfacn,wsfacn=NULL,sfac="s",SStype=3,spss=F,dvnam=NA) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  if (is.na(dvnam)) dvnam <- names(dat)[dim(dat)[2]]
  dat <- dat[,c(sfac,bsfacn,wsfacn,dvnam)]
  if (length(wsfacn)>0) dat <- dat[do.call(order,dat[,c(sfac,wsfacn)]),]
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  snams=levels(dat[,sfac]); ns=length(snams)
  nifac=length(wsfacn)
  lev1s=unlist(lapply(lapply(dat[,wsfacn],levels),function(x){x[1]}))
  bsfacs=dat[apply(dat[,wsfacn,drop=F],1,function(x){all(x==lev1s)}),bsfacn,drop=F]  
  if ( nifac>0 ) {
    idata=data.frame(dat[dat[,sfac]==snams[1],wsfacn])
    names(idata)=wsfacn
    for (i in wsfacn) 
      if (i==wsfacn[1]) ifacn=as.character(idata[,1]) else 
        ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
    e.mv=matrix(unlist(
      tapply(dat[,dvnam],dat[,wsfacn[length(wsfacn):1]],function(x){x})),
                ncol=length(ifacn),dimnames=list(snams,ifacn))
    print(summary(Anova(
      lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),bsfacs),
      idata=idata,type=SStype,
      idesign=formula(paste("~",paste(wsfacn,collapse="*")))),multivariate=F))
  } else {
    e.mv <- cbind(y=dat[,dvnam]) 
    print(Anova(lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),
                   bsfacs),type=3))
  }
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),bsfacs,e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}

# Error Bars for Repeated Measures ----------------------------------------

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}



# Get Normalized Standard Errors for One Subject --------------------------

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Normalize Data Within Groups --------------------------------------------

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}


# Cheater error bars ------------------------------------------------------

SE <- function(dependent_variable) {
  sd(dependent_variable)/sqrt(length(dependent_variable))
}

lower <- function(dependent_variable) {
  mean(dependent_variable)-2*SE(dependent_variable)
}
upper <- function(dependent_variable) {
  mean(dependent_variable)+2*SE(dependent_variable)
}




# scol = subject column (numeric or name, same for following)
# dvcol = dependent variable column
# ws and bs name corresponding factors in dat
# mixedAnova=function(dat,ws,bs=NA,scol=1,dvcol=dim(dat)[2],
#   SStype=3,icontrasts=c("contr.sum", "contr.poly"),
#   savedf=F) {
#   dat=dat[do.call(order,dat[,ws,drop=F]),]
#   snams=levels(dat[,scol]); ns=length(snams)
#   dvnam=names(dat)[dvcol]
#   nifac=length(ws)
#   idata=data.frame(dat[dat[,scol]==snams[1],ws])
#   names(idata)=ws  
#   for (i in ws) 
#     if (i==ws[1]) ifacn=as.character(idata[,1]) else 
#                   ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
#   facnr=ws[nifac:1]
#   e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
#               ncol=length(ifacn),dimnames=list(snams,ifacn))
#   if (length(bs)==1) bsf=bs else bsf=paste(bs,collapse="*")
#   if (any(is.na(bs))) {
#     form=formulua("e.mv~1")
#     tmp <- e.mv
#   } else {
#     form=formula(paste("e.mv ~",bsf))
#     trans=snams
#     names(trans)=snams
#     for (i in bs) {
#       for (j in snams) trans[j] <- as.character(dat[dat[,scol]==j,i][1])
#       if (i==bs[1]) bsfac <- factor(trans[dimnames(e.mv)[[1]]]) else
#         bsfac <- cbind.data.frame(bsfac,factor(trans[dimnames(e.mv)[[1]]]))
#     }
#     tmp <- cbind.data.frame(bsfac,e.mv)
#     names(tmp)[1:length(bs)] <- bs
#   }
#   summary(Anova(lm(form,tmp),
#     idata=idata,type=SStype,icontrasts=icontrasts,
#     idesign=formula(paste("~",paste(ws,collapse="*"))) ),
#     multivariate=FALSE)
#   invisible(tmp)
# }
