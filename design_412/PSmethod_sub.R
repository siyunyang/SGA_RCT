#' Fitting propensity scores with different models
#'
#' The function \code{PSmethod} is an internal function to estimate the propensity scores given a specified model through formula.
#' It is bulit into function \code{Sumstat}, \code{PStrim} and \code{PSweight}.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param method a character to specify the method for propensity model. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param data an optional data frame containing the variables in the propensity score model.
#' @param ncate a numeric to specify the number of treatment groups present in the given data.
#' @param ... further arguments passed to or from other methods.
#'
#' @details  A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable and \code{terms} is a series of terms which specifies a linear predictor for \code{treatment}. \code{ps.formula} by default specifies generalized
#' linear models given the default argument \code{method = "glm"}.  It fits the logistic regression when \code{ncate = 2},and multinomial
#' logistic regression when \code{ncate > 2}. The argument \code{method} allows user to choose
#' model other than glm to fit the propensity score models. We have included \code{LASSO}, \code{pLASSO}, \code{gbm} and \code{SuperLearner}.
#' Additional argument in them can be supplied through the \code{...} argument. Note that SuperLearner does not handle multiple groups and the current version of multinomial
#' logistic regression is broken in gbm. We suggest user to use them with extra caution. Please refer to the user manual of the \code{gbm} and \code{SuperLearner} packages for all the
#' allowed arguments.
#'
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ e.h}}{a data frame of estimated propensity scores.}
#'
#' }
#'
#' @export
#'
#' @examples
#' # the propensity model
#' ps.formula <- trt~cov1+cov2+cov3+cov4+cov5+cov6
#' psfit <- PSmethod(ps.formula = ps.formula,data = psdata,ncate=3)
#'
#'
#'
#'ps.formula specify all covariates, may include confounders, subgroup variables and their interactions
#'Connect plot: default ALL = FALSE (covariate - subgroup -interaction); All= TRUE (main effect including confounders and subgroups)
#'remove trimming
#'
#'
#'
#'

library(rms)
library(gbm)
library(glmnet)
library(ranger)
PSmethod_sub <-function(ps.formula=ps.formula, subgroup=NULL, method="glm", data=data, ...){

  ps.formula<-as.formula(ps.formula)
  zname<-all.vars(ps.formula)[1]
  facz<-as.factor(data[,zname])
  #creat a dictionary for the original and recoded values in Z
  dic<-levels(facz)
  
  
  if (method=="glm"){
  ############## logistic #############################################################
    #change z to 0/1
    dataps<-data
    dataps[,zname]<- as.numeric(facz)-1

    fitglm <- glm(formula = ps.formula, data=dataps,family = binomial(link = "logit"))
    e.h <- fitglm$fitted.values
    e.h <- cbind(1-e.h,e.h)
  }else if (method=="gbm"){
  ############## gbm ###################################################################
    
    if (exists("distribution")){
      if (!distribution %in% c("bernoulli","adaboost","multinomial")){
        stop("only bernoulli, adaboost, or multinomial distributions in 'gbm' are supported in propensity score models of PSweight")
      }
    }

    if (exists("var.monotone")) var.monotone=var.monotone else var.monotone=NULL
    if (exists("weights")) weights=weights else weights=NULL
    if (exists("n.trees")) n.trees=n.trees else n.trees=100
    if (exists("interaction.depth")) interaction.depth=interaction.depth else interaction.depth=1
    if (exists("n.minobsinnode")) n.minobsinnode=n.minobsinnode else n.minobsinnode=10
    if (exists("shrinkage")) shrinkage=shrinkage else shrinkage=0.1
    if (exists("bag.fraction")) bag.fraction=bag.fraction else bag.fraction=0.5
    if (exists("train.fraction")) train.fraction=train.fraction else train.fraction=1
    if (exists("cv.folds")) cv.folds=cv.folds else cv.folds=0
    if (exists("class.stratify.cv ")) class.stratify.cv =class.stratify.cv  else class.stratify.cv=NULL
    if (exists("n.cores ")) n.cores =n.cores  else n.cores=NULL
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")

 
    #change z to 0/1
    dataps<-data
    dataps[,zname]<- as.numeric(facz)-1

    if (exists("distribution")) {
      if (!distribution %in% c("adaboost","bernoulli")) {
        distribution<-"bernoulli"
        warning("supplied unsupported distribution for binary outcome in gbm; reset to bernoulli")
      }
    }else{
      distribution<-"bernoulli"
    }

    fitgbm <- gbm::gbm(formula = ps.formula, data=dataps,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                     interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                     train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                     class.stratify.cv = class.stratify.cv, n.cores = n.cores)

    e.h<-exp(fitgbm$fit)/(1+exp(fitgbm$fit))
    e.h<-cbind(1-e.h,e.h)
  }else if (method =='rcs'){
  ############## fcs #############################################################
    #change z to 0/1
    dataps<-data
    dataps[,zname]<- as.numeric(facz)-1
    
    fi_rcs <- lrm(formula = ps.formula, data=data, x=TRUE, y =TRUE) 
    e.h <- predict(fi_rcs, type="fitted.ind")
    e.h <- cbind(1-e.h,e.h)
  }else if(method =='pLASSO'){
    ############## pLASSO #############################################################
    #change z to 0/1
    dataps<-data
    z<- as.numeric(facz)-1
    dataps[,zname]<-z
    
    fullmatrix <- model.matrix(ps.formula, data)
    
    # only penalize interactions
    # create an idx for penalty.factor 
    idx <- rep(0, dim(fullmatrix)[2])
    idx[grep('\\:',colnames(fullmatrix))] <- 1
    
    fitLASSO <- cv.glmnet(y=factor(z), x=fullmatrix, penalty.factor=idx, family="binomial", maxit=50000)
    nonzero_coef <- rownames(coef(fitLASSO, s='lambda.min'))[which(coef(fitLASSO, s='lambda.min')!=0)][-1]
    if(length(nonzero_coef)==0){
      warning("no coefficient ")
      fitpLASSO <- glm(factor(z)~1, family = binomial(link = "logit"))
    }else{ 
      fitpLASSO <- glm(factor(z)~., data=data.frame(fullmatrix[,nonzero_coef]), family = binomial(link = "logit"))
    }
    e.h <- fitpLASSO$fitted.values
    e.h <- cbind(1-e.h,e.h)
  }else if (method=="RFs"){
    ############## RFs ###################################################################
    
    if (exists("num.trees")) num.trees=num.trees else num.trees=500
    if (exists("num.threads")) num.threads=num.threads else num.threads=10
    if (exists("seed")) seed=seed else seed=123
    if (exists("verbose")) warning("verbose argument set to F for RFs in PSweight")
    if (exists("probability")) warning("probability argument set to T for RFs in PSweight")
    
    
    #change z to 0/1
    dataps<- data
    dataps[,zname]<- as.factor(dataps[,zname])
    
    fitrfs <- ranger::ranger(formula = ps.formula, data = dataps, num.trees = num.trees, replace = T,  probability = T, verbose = F,
                     num.threads =num.threads, importance = "none", seed = seed)
    # outbag prediction
    e.h <- fitrfs$predictions[,2]
    e.h <- cbind(1-e.h, e.h)
  }
  
  colnames(e.h)<-dic

  return(e.h=e.h)
}





















