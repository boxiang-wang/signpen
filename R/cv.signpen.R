#' Cross-validation for signpen
#' 
#' Does k-fold cross-validation for signpen and returns a
#' value for \code{lambda}. This function is modified based on the \code{cv}
#' function from the \code{glmnet} package.
#' 
#' The function runs \code{\link{signpen}} \code{nfolds}+1 times; the first to
#' get the \code{lambda} sequence, and then the remainder to compute the fit
#' with each of the folds omitted. The average error and standard deviation
#' over the folds are computed.
#' 
#' @aliases cv.signpen cv.signlog cv.signls
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#'   is an observation vector.
#' @param y response variable. This argument should be quantitative for least squares, 
#'   and a two-level factor for logistic regression.
#' @param lambda optional user-supplied lambda sequence; default is
#'   \code{NULL}, and \code{\link{signpen}} chooses its own sequence.
#' @param pred.loss loss to use for cross-validation error. Valid options are:
#'   \itemize{ \item \code{"loss"} for classification, margin based loss
#'   function.  \item \code{"misclass"} for classification, it gives
#'   misclassification error.  \item \code{"L1"} for regression, mean square
#'   error used by least squares regression \code{loss="ls"}, it measure the
#'   deviation from the fitted mean to the response.  \item \code{"L2"} for
#'   regression, mean absolute error used by least squares regression
#'   \code{loss="ls"}, it measure the deviation from the fitted mean to the
#'   response.  } Default is \code{"loss"}.
#' @param nfolds number of folds; default is 5. 
#'   Smallest value allowable is \code{nfolds=3}.
#' @param foldid an optional vector of values between 1 and \code{nfold}
#'   identifying what fold each observation is in. 
#'   If supplied, \code{nfold} can be missing.
#' @param \dots other arguments that can be passed to signpen.
#' @return an object of class \code{\link{cv.signpen}} is returned, which is a
#'   list with the ingredients of the cross-validation fit.  \item{lambda}{the
#'   values of \code{lambda} used in the fits.} \item{cvm}{the mean
#'   cross-validated error - a vector of length \code{length(lambda)}.}
#'   \item{cvsd}{estimate of standard error of \code{cvm}.} \item{cvupper}{upper
#'   curve = \code{cvm+cvsd}.} \item{cvlower}{lower curve = \code{cvm-cvsd}.}
#'   \item{name}{a text string indicating type of measure (for plotting
#'   purposes).} \item{signpen.fit}{a fitted \code{\link{signpen}} object for the
#'   full data.} \item{lambda.min}{The optimal value of \code{lambda} that gives
#'   minimum cross validation error \code{cvm}.} \item{lambda.1se}{The largest
#'   value of \code{lambda} such that error is within 1 standard error of the
#'   minimum.}
#' @keywords perceptron penalty sign consistency
#' @export
cv.signpen = function(x, y, lambda=NULL, 
  pred.loss=c("misclass", "loss"), nfolds=5, foldid, ...) {
  ####################################################################
  ## data setup
  y = as.numeric(drop(y))
  x = as.matrix(x)
  x.row = as.integer(NROW(x))
  if (length(y) != x.row) 
    stop("x and y have different number of observations.")  
  ####################################################################
  signpen.object = signpen(x, y, lambda=lambda, ...)
  lambda = signpen.object$lambda
  if (missing(foldid)) 
    foldid = sample(rep(seq(nfolds), 
      length=x.row)) else nfolds = max(foldid)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended.")
  outlist = as.list(seq(nfolds))
  ## fit the model nfold times and save them
  for (i in seq(nfolds)) {
    which = foldid == i
    outlist[[i]] = signpen(x=x[!which, , drop=FALSE], 
      y=y[!which], lambda=lambda, ...)
  }
  if (missing(pred.loss)) 
    pred.loss = "default" else pred.loss = match.arg(pred.loss)
  ## select the lambda according to predmat
  fun = paste("cv", class(signpen.object)[[2]], sep=".")
  cvstuff = do.call(fun, 
    list(outlist, lambda, x, y, foldid, pred.loss))
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  cvname = cvstuff$name
  out = list(lambda=lambda, cvm=cvm, cvsd=cvsd, 
    cvupper=cvm+cvsd, cvlower=cvm - cvsd, 
    name=cvname, signpen.fit=signpen.object)
  obj = c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) = "cv.signpen"
  obj
} 

cv.signlog = function(outlist, lambda, x, y, foldid, pred.loss) {

  if (pred.loss == "default") pred.loss = "misclass" 
  if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
    warning("Only 'misclass' and 'loss' are available for logit; 'misclass' used.")
    pred.loss = "misclass"
  }
  typenames = c("misclass" = "mis-classification error",
                "loss"     = "binomial deviance")
  ### Turn y into c(-1,1)
  y = c(-1, 1)[as.numeric(y)]
  nfolds = max(foldid)
  predmat = matrix(NA, length(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, x[which, , drop=FALSE], type="link")
    nlami = length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  cvraw = switch(pred.loss, 
    "loss"     = logitloss(y * predmat),
    "misclass" = (y != ifelse(predmat > 0, 1, -1)))
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw
    cvn = cvob$N
  } else cvn = length(y) - colSums(is.na(predmat))    
  cvm = colMeans(cvraw, na.rm=TRUE)
  cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, 
    na.rm=TRUE)/(cvn - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 

cv.signls = function(outlist, lambda, x, y, foldid, pred.loss) {
  if (pred.loss == "default") pred.loss = "loss" 
  if (!match(pred.loss, c("loss"), FALSE)) {
    warning("Only 'loss' available for least squares regression; 'loss' used")
    pred.loss = "loss"
  }
  typenames = c(loss = "least-squared loss")
  y = as.numeric(y)
  nfolds = max(foldid)
  predmat = matrix(NA, length(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    preds = predict(fitobj, x[which, , drop=FALSE], type="link")
    nlami = length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] = preds
    nlams[i] = nlami
  }
  cvraw = (y - predmat)^2
  if (length(y)/nfolds >= 3) {
    cvob = cvcompute(cvraw, foldid, nlams)
    cvraw = cvob$cvraw
    cvn = cvob$N
  } else cvn = length(y) - colSums(is.na(predmat))    
  cvm = colMeans(cvraw, na.rm=TRUE)
  cvsd = sqrt(colMeans(scale(cvraw, cvm, FALSE)^2, 
    na.rm=TRUE)/(cvn - 1))
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
