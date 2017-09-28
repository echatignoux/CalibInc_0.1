##' Prediction of incidence from proxy data using a calibration approach
##'
##' \code{CalibInc} functions takes a calibration model and a table with proxy data as
##' arguments and returns a table with predictions and standard errors
##'
##' @param mod Calibration model (type H/I ratio), fitted with
##'   \code{glmer} or \code{glmmPQL}
##' @param pred Data frame containing proxy data to be calibrated. If \code{NULL}, try to use the
##' data.frame used for model evaluation (searched in .GlobalEnv).
##' @param weight Weight vector (size 1 for an uniform weight, same size as
##'   \code{nrow(pred)} otherwise). \code{weight} can be given as the name of a
##'   column in \code{pred}.
##' @param aggregate Aggregate the predictions according to a combination of variables given as a formula
##'    (see examples). If \code{NULL} aggregation is done according to the combinations of all co-variables
##'    specified in the model formula.
##' @param keep.Vp Set to \code{TRUE} to keep variances covariance matrix of the predictions
##' @return A data frame containing predictions values (column \code{pred}) and se (column \code{se})
##'    at the levels defined in \code{aggregate}. If keep.Vp=\code{TRUE}, the
##'    variances covariance matrix of the predictions is also returned as an attribute.
##' @author Edouard Chatignoux
##' @export
##' @examples
##' library(lme4)
##' library(splines)
##' library(tidyverse)
##' data(lopm.CalibSet)
##' data(lopm.Fr)
##' ## Run the calibration model for lop cancer incidence
##' ## and hospitalisation data in men, 2007-2011
##' k<-Hmisc::wtd.quantile(lopm.CalibSet$age,weights=lopm.CalibSet$C,p=0.5)%>%as.numeric
##' form.calib<-substitute(H~offset(log(C))+
##'                        ns(age,knots=k, Boundary.knots = range(age)+c(5,-5))+
##'                        (1|dist),
##'                        list(k=k))%>%as.formula
##' mod.calib<-glmer(form.calib,data=lopm.CalibSet%>%filter(C>0),family="poisson", nAGQ =20)
##' ## Predict the total number of incident cases by district
##' CalibInc(mod.calib,pred=lopm.Fr,aggregate=~dist)
##' ## Predict the total number of incident cases by age and district
##' CalibInc(mod.calib,pred=lopm.Fr,aggregate=~dist+age)
CalibInc<-  function (mod, pred = NULL,
                      weight = 1,
                      aggregate=NULL,
                      keep.Vp = FALSE)  {

  vars <- all.vars(formula(mod))
  nb.can_ <- vars[attr(terms(mod),"offset")]
  nb.bma_ <- vars[1]

  ag<-aggregate
  if (is.null(ag))
    ag<-paste("~",paste(rev(setdiff(vars,c(nb.bma_,nb.can_))),collapse=":"))%>%as.formula
  else if(ag!=~1) ag<-paste("~",paste(all.vars(ag),collapse=":"))%>%as.formula
  if(ag!=~1) ag<-update(ag,~.-1)

  if (is.null(pred))
    if (class(mod)[1] == "glmerMod")
      pred <- try(eval(attr(mod, "call")$data,.GlobalEnv), silent = TRUE)
  else pred <- try(eval(mod$call$data, .GlobalEnv), silent = TRUE)
  if (class(pred)[1]=="try-error")
    stop("The GLMM data set has changed since evaluation. Specify a prediction dataset trough the pred argument.")
  pred<-pred%>%dplyr::tbl_df()%>%dplyr::ungroup()
  if (any(is.na(match(c(all.vars(ag),setdiff(all.vars(nobars(formula(mod))),nb.can_)),names(pred)))))
    stop("Variables specified in aggregate or in the model are missing from the pred table.")
  if(ag!=~1)
    pred<-pred%>%dplyr::arrange_(all.vars(ag))

  if (!(class(mod)[1] %in% c("glmerMod","glmmPQL")))
    stop("Only glmerMod or glmmPQL object are allowed in mod.")

  if (is.name(match.call()$weight)) w<-pred[,deparse(substitute(weight))]%>%unlist
  else if (is.character(weight)) w<-pred[,weight]%>%unlist
  else w<-weight


  W <- matrix(w, nrow = 1, ncol = nrow(pred))
  if (class(mod)[1] == "glmerMod") {
    beta <- fixef(mod)
    V <- vcov(mod)
    sigma <- sigma(mod)
    sb <- mod@theta^2
  }
  if (class(mod)[1] == "glmmPQL") {
    beta <- mod$coefficients$fixed
    V <- mod$varFix
    sigma <- mod$sigma^2
    sb <- mod$modelStruct$reStruct[[1]][1] * sigma
  }
  pred[, nb.can_] <- 1
  Xp <- model.matrix(nobars(formula(mod)), pred)
  alpha <- c(beta, sb)
  Z <- cbind(Xp, Xp[, 1]/2)
  Sa <- bdiag(list(V, 0))
  S <- Z %*% Sa %*% t(Z)
  eZa <- exp(-Z %*% alpha)
  BMA <- dplyr::collect(dplyr::select(pred, dplyr::contains(nb.bma_)))[[1]]
  P <- as.matrix(BMA * eZa, ncol = 1)
  Ds <- matrix(diag(as.matrix(S)))
  Vp <- Matrix::Diagonal(x = sigma * P * (exp(-Z %*% alpha +
                                        2 * Ds) + 1)) + (exp(S + sb) - 1) * exp(Ds/2) %*%
    t(exp(Ds/2)) * (Matrix::Diagonal(x = P) + P %*% t(P))

  nbbma <- pred[,nb.bma_][[1]]%>% as.vector

  f <- function(x) factor(x, levels = unique(x))
  Xag<-pred%>%dplyr::mutate_at(all.vars(ag),funs(f))
  Xag <- model.matrix(ag, data = droplevels(Xag))
  Xag <- Xag[,colSums(Xag)!=0]
  Xag <- t(Xag) %*% Matrix::Diagonal(x = as.vector(W))
  if(ag!=~1)
    dtPred = pred%>%dplyr::select_(.dots =all.vars(ag))%>%unique%>%
      dplyr::mutate(nbbma = Xag %*% nbbma %>% as.vector,
                    pred = Xag %*% P %>% as.vector,
                    se = sqrt(Matrix::diag(Xag %*%Vp %*% Matrix::t(Xag))))%>%
      dplyr::rename_(.dots=setNames(list("nbbma"), list(nb.bma_)))
  else
    dtPred = pred%>%dplyr::select_(.dots =all.vars(ag)) %>%
      dplyr::mutate(nbbma = Xag %*% nbbma %>% as.vector,
                    pred = Xag %*% P %>% as.vector,
                    se = sqrt(Matrix::diag(Xag %*%Vp %*% Matrix::t(Xag))))%>%
      dplyr::rename_(.dots=setNames(list("nbbma"), list(nb.bma_)))%>%
      unique

  attr(dtPred, "Sb") <- sb * sigma
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


