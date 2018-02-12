##' Predict fixed effects part of a gam model with random-effect
##'
##'
##' @param mod A model fitted with gam
##' @param newdata Prediction data set
##' @param type \code{lpmatrix} to get the model matrix, \code{response} for the response
##' @return A matrix or a vector
##' @author Edouard Chatignoux
##' @export
##' @keywords internal
##' @examples
predict_gam_fixed<-function(mod,newdata=NULL,type="lpmatrix"){
  if (is.null(newdata)) newdata<-mod$model
  st<-sapply(mod$smooth,function(m)m$term)
  sre<-sapply(mod$smooth,function(m)ifelse(is.null(m$random),FALSE,m$random))
  ret<-st[sre]
  wre<-lapply(mod$smooth[sre],function(sre)seq(sre$first.para,sre$last.para))%>%do.call("c",.)
  Xs<-lapply(1:length(sre),function(i)
    if (!sre[i]) mgcv::PredictMat(mod$smooth[[i]], newdata))%>%
    do.call("cbind",.)
  drop_sm<-which((attr(terms(mod),"variables")%>%as.character)[-c(1,2,attr(terms(mod), "offset")+1)] %in% st)
  drop_sm<-paste0(".~.-",paste0(attr(terms(formula(mod)),"term.labels")[drop_sm],collapse="-"))

  pred<-newdata
  nst<-setdiff(all.vars(formula(mod)),st)
  for (v in setdiff(nst,names(pred)))  pred[,v]<-1
  Xf<-model.matrix(update(formula(mod),drop_sm),pred%>%data.frame)
  Xp<-cbind(Xf,Xs)

  if (type=="lpmatrix") return(Xp)
  if (type=="response"){
    beta <- mod$coefficients[-wre]
    V <- mod$Vp[-wre,-wre]
    pred<-Xp%*%beta%>%as.vector
    se.pred<-sqrt(diag(Xp%*%V%*%t(Xp)))
    newdata%>%tbl_df%>%
      mutate(pred=exp(pred),
             low=pred*exp(-1.96*se.pred),
             up=pred*exp(+1.96*se.pred))
    }
}


##' Prediction of incidence from proxy data using a calibration approach
##'
##' \code{CalibInc} functions takes a calibration model and a table with proxy data as
##' arguments and returns a table with predictions and standard errors
##'
##' @param mod Calibration model (type H/I ratio), fitted with
##'   \code{glmer}, \code{glmmPQL}, \code{gam} or \code{gamm}
##' @param pred Data frame containing proxy data to be calibrated. If
##'   \code{NULL}, try to use the data.frame used for model evaluation
##'   (searched in .GlobalEnv).
##' @param weight Weight vector (size 1 for an uniform weight, same
##'   size as \code{nrow(pred)} otherwise). \code{weight} can be given
##'   as the name of a column in \code{pred}.
##' @param aggregate Aggregate the predictions according to a
##'   combination of variables given as a formula (see examples). If
##'   \code{NULL} aggregation is done according to the combinations of
##'   all co-variables specified in the model formula.
##' @param keep.Vp Set to \code{TRUE} to keep variances covariance
##'   matrix of the predictions
##' @return A data frame containing predictions values (column
##'   \code{pred}) and se (column \code{se}) at the levels defined in
##'   \code{aggregate}. If keep.Vp=\code{TRUE}, the variances
##'   covariance matrix of the predictions is also returned as an
##'   attribute.
##' @author Edouard Chatignoux
##' @export
##' @examples
##' library(lme4)
##' library(splines)
##' library(tidyverse)
##' data(lopm.CalibSet)
##' data(lopm.alldist)
##' ## Run the calibration model for lop cancer incidence
##' ## and hospitalisation data in men, 2007-2011
##' k<-Hmisc::wtd.quantile(lopm.CalibSet$age,weights=lopm.CalibSet$C,p=0.5)%>%as.numeric
##' form.calib<-substitute(H~offset(log(C))+
##'                        ns(age,knots=k, Boundary.knots = range(age)+c(5,-5))+
##'                        (1|dist),
##'                        list(k=k))%>%as.formula
##' mod.calib<-glmer(form.calib,data=lopm.CalibSet%>%filter(C>0),family="poisson", nAGQ =20)
##' ## Predict the total number of incident cases by district
##' CalibInc(mod.calib,pred=lopm.alldist,aggregate=~dist)
##' ## Predict the total number of incident cases by age and district
##' CalibInc(mod.calib,pred=lopm.alldist,aggregate=~dist+age)
CalibInc<-function (mod, pred = NULL, weight = 1, aggregate = NULL, keep.Vp = FALSE)
{
  f <- function(x) factor(x, levels = unique(x))
  if (!(class(mod)[1] %in% c("glmerMod", "glmmPQL","gam","gamm")))
    stop("Only gam(m), glmerMod and glmmPQL object are allowed in mod.")
  ## Model parameters
  mod_<-NULL
  if (class(mod)[1] == "glmerMod") {
    beta <- lme4::fixef(mod)
    V <- stats::vcov(mod)
    sigma <- sigma(mod)
    sb <- mod@theta^2
    re<-names(summary(mod)$ngrps)
  }
  if (class(mod)[1] == "glmmPQL") {
    beta <- mod$coefficients$fixed
    V <- mod$varFix
    sigma <- mod$sigma^2
    sb <- mod$modelStruct$reStruct[[1]][1] * sigma
    re<-names(mod$groups)
  }
  if (class(mod)[1] == "gam") {
    sre<-sapply(mod$smooth,function(m)ifelse(is.null(m$random),FALSE,m$random))
    re<-sapply(mod$smooth,function(m)m$term)[sre]
    wre<-grep(re,names(mod$coefficients))
    beta <- mod$coefficients[-wre]
    V <- vcov(mod)[-wre,-wre]
    sigma <- summary(mod)$dispersion
    tmp<-tempfile()
    sink(tmp)
    if (mod$method=="REML")
      sb <- mgcv::gam.vcomp(mod)[sre,1]^2
    else sb <- mgcv::gam.vcomp(mod)[sre]^2
    sink()
    file.remove(tmp)
  }
  if (class(mod)[1] == "gamm") {
    beta <- mod$gam$coefficients
    V <- vcov(mod$gam)
    sigma <- summary(mod$gam)$dispersion
    sb <- mod$lme$modelStruct$reStruct[[1]][1] * sigma
    re<-names(mod$lme$groups)[-1]
    mod_<-mod
    mod<-mod$gam
  }
  if (length(re)>1)
    stop("Only one random-effect is allowed.")

  ## Covariates and aggregation variables
  vars <- all.vars(formula(mod))
  nb.can_ <- vars[attr(terms(mod), "offset")]
  if (class(mod)[1] == "gam" & !is.null(mod_))
    if (!is.null(attr(terms(mod), "offset")))
      nb.can_ <- as.character(attr(terms(mod),"variables")[[attr(terms(mod), "offset")+1]])[2]
  nb.bma_ <- vars[1]
  covars <- as.character(sapply(attr(terms(formula(mod)), "term.labels"),
                   function(tl) all.vars(as.formula(paste0("~", tl)))[1]))
  if (class(mod)[1] == "glmmPQL")
    covars <- c(covars, names(mod$groups))
  if (class(mod)[1] == "gam")
    covars <- c(covars, names(mod_$lme$groups)[-1])
  ag <- aggregate
  if (is.null(ag))
    ag <- stats::as.formula(paste("~", paste(covars, collapse = ":")))
  else if (ag != ~1)
    ag <- stats::as.formula(paste("~", paste(all.vars(ag), collapse = ":")))
  if (ag != ~1)
    ag <- update(ag, ~. - 1)
  ## Prediction table
  if (is.null(pred)){
    if (class(mod)[1] == "glmerMod")
      pred <- try(base::eval(attr(mod, "call")$data, .GlobalEnv),
                  silent = TRUE)
    else if ((class(mod)[1] == "gam" & is.null(mod_)))
      pred<-mod$model
    else if ((class(mod)[1] == "gam" & !is.null(mod_)))
      pred<-data.frame(mod$model%>%dplyr::select(-g:-X))
    else pred <- try(base::eval(mod$call$data, .GlobalEnv), silent = TRUE)
    }
  if (class(pred)[1] == "try-error")
    stop("The GLMM data set has changed since evaluation. Specify a prediction dataset trough the pred argument.")
  pred <- pred %>% dplyr::tbl_df() %>% dplyr::ungroup()
  if (ag != ~1)
    pred <- pred %>% dplyr::arrange_(rev(all.vars(ag)))
  pred[, nb.can_] <- 1
  ## Weight
  if (is.name(match.call()$weight))
    w <- pred[, deparse(substitute(weight))] %>% unlist
  else if (is.character(weight))
    w <- pred[, weight] %>% unlist
  else w <- weight
  W <- Matrix::Diagonal(x = as.vector(matrix(w, nrow = 1, ncol = nrow(pred))))
  if (any(is.na(match(re, names(pred)))))
    stop("Prediction table must contain the district variable.")
  if (any(is.na(match(c(all.vars(ag), covars, nb.bma_), names(pred)))))
    stop("Variables specified in aggregate or in the model are missing from the pred table.")
  ## Prediction and aggregation matrix
  if (class(mod)[1] != "gam") Xp <- model.matrix(lme4::nobars(formula(mod)), pred)
  else if (!is.null(mod_)) Xp <- predict(mod, newdata=pred,type="lpmatrix")
  else if (class(mod)[1] == "gam") Xp <- predict_gam_fixed(mod, newdata=pred)
  if (length(unique(pred[,re][[1]]))>1)
    Xdep <- model.matrix(as.formula(paste0("~-1+", re)),
                         data = pred %>%
                           dplyr::mutate_at(re, funs(f))) %>% Matrix::Matrix(.,sparse = T)
  else   Xdep <- model.matrix(as.formula(~1),
                              data = pred) %>% Matrix::Matrix(.,sparse = T)
  Xag <- pred %>% dplyr::mutate_at(all.vars(ag), funs(f))
  Xag <- model.matrix(ag, data = droplevels(Xag))
  Xag <- Xag[, colSums(Xag) != 0]
  Xag <- Matrix::t(Xag) %*% W
  ## Predictions and prediction-variance
  alpha <- c(beta, sb)
  Z <- cbind(Xp, Xp[, 1]/2)%>%Matrix::Matrix(.,sparse=T)
  Sa <- Matrix::bdiag(list(V, 0))
  S <- Z %*% Sa %*% Matrix::t(Z)
  eZa <- exp(-Z %*% alpha)
  BMA <- dplyr::collect(dplyr::select(pred, dplyr::contains(nb.bma_)))[[1]]
  P <- BMA * eZa
  Ds <- Matrix::Matrix(diag(as.matrix(S)))
  eDs<-Matrix::tcrossprod(exp(Ds/2),exp(Ds/2))
  PtP<-Matrix::tcrossprod(P,P)
  XtX<-Matrix::tcrossprod(Xdep,Xdep)
  Vp <- Matrix::Diagonal(x = (sigma * P * (exp(-Z %*% alpha +2 * Ds) + 1))[,1]) +
    (exp(S) - 1) * eDs * (Matrix::Diagonal(x = P[,1]) + PtP)+
    (exp(S)*(exp(sb) - 1))*eDs*(PtP)*(XtX)
  ## Output
  if (ag != ~1)
    dtPred = pred %>% dplyr::select_(.dots = all.vars(ag)) %>%
      unique %>% dplyr::mutate(nbbma = Xag %*% BMA %>%as.vector,
                               pred = Xag %*% P %>% as.vector,
                               se = sqrt(Matrix::diag(Xag %*%Vp %*% Matrix::t(Xag)))) %>%
      dplyr::rename_(.dots = setNames(list("nbbma"),list(nb.bma_)))
  else dtPred = pred %>% dplyr::select_(.dots = all.vars(ag)) %>%
         dplyr::mutate(nbbma = Xag %*% BMA %>% as.vector,
                       pred = Xag %*%P %>% as.vector,
                       se = sqrt(Matrix::diag(Xag %*% Vp %*%Matrix::t(Xag)))) %>%
         dplyr::rename_(.dots = setNames(list("nbbma"),list(nb.bma_))) %>% unique
  attr(dtPred, "Sb") <- sqrt(sb)
  if (keep.Vp)
    attr(dtPred, "Vp") <- Vp
  dtPred
}


##' Confidence interval for log-normal distributed observations
##'
##' \code{LogNormPI} calcultes confidence interval for log-normal distributed observations
##'
##' @param data A data frame with log-normal observations
##' @param pred Log-normal observation
##' @param se Standard error of \code{pred}
##' @param level Confidence interval level
##' @return
##' @author Edouard Chatignoux
##' @export
##' @examples
##' dt<-data_frame(m=2,sd=1)
##' LogNormPI(data=dt,m,sd)
LogNormPI<-function(data,pred=pred,se=se,level=.95){
  z<-qnorm(1-(1-level)/2)
  data$pred_<-data[deparse(substitute(pred))]%>%unlist
  data$se_<-data[deparse(substitute(se))]%>%unlist
  data%>%dplyr::mutate(cv=se_/pred_,
                low=pred_*sqrt(cv^2+1)*exp(-z*sqrt(log(cv^2+1))),
                up=pred_*sqrt(cv^2+1)*exp(+z*sqrt(log(cv^2+1))))%>%
    dplyr::select(-cv,-pred_,-se_)}






