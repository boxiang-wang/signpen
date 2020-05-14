#' get coefficients or make coefficient predictions from a "cv.signpen" object.
#' 
#' This function gets coefficients or makes coefficient predictions from a
#'   cross-validated \code{signpen} model, using the stored \code{"signpen.fit"}
#'   object, and the optimal value chosen for \code{lambda}.
#' 
#' This function makes it easier to use the results of cross-validation to get
#'   coefficients or make coefficient predictions.
#' 
#' @param object fitted \code{\link{cv.signpen}} object.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#'   predictions are required. Default is the value \code{s="lambda.1se"} stored
#'   on the CV \code{object}, it is the largest value of \code{lambda} such that
#'   error is within 1 standard error of the minimum. Alternatively
#'   \code{s="lambda.min"} can be used, it is the optimal value of \code{lambda}
#'   that gives minimum cross validation error \code{cvm}. If \code{s} is
#'   numeric, it is taken as the value(s) of \code{lambda} to be used.
#' @param \dots not used. Other arguments to predict.
#' @return The coefficients at the requested values for \code{lambda}.
#' @keywords perceptron penalty sign consistency
#' @export
coef.cv.signpen = function(object, s=c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) 
    lambda = s else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
  } else stop("Invalid form for s")
  coef(object$signpen.fit, s=lambda, ...)
}

#' get coefficients or make coefficient predictions from an "signpen" object.
#' 
#' Computes the coefficients at the requested values for \code{lambda} from a
#'   fitted \code{\link{signpen}} object.
#' 
#' \code{s} is the new vector at which predictions are requested. If \code{s}
#'   is not in the lambda sequence used for fitting the model, the \code{coef}
#'   function will use linear interpolation to make predictions. The new values
#'   are interpolated using a fraction of coefficients from both left and right
#'   \code{lambda} indices.
#' 
#' @param object fitted \code{\link{signpen}} model object.
#' @param s value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the entire sequence used to create the model.
#' @param \dots not used. Other arguments to predict.
#' @return The coefficients at the requested values for \code{lambda}.
#' @keywords perceptron penalty sign consistency
#' @export
coef.signpen = function(object, s=NULL, ...) {
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    if(length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
        nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  return(nbeta)
}

#' make predictions from a "signpen" object.
#' 
#' Similar to other predict methods, this functions predicts fitted values and
#'   class labels from a fitted \code{\link{signpen}} object.
#' 
#' \code{s} is the new vector at which predictions are requested. If \code{s}
#'   is not in the lambda sequence used for fitting the model, the \code{predict}
#'   function will use linear interpolation to make predictions. The new values
#'   are interpolated using a fraction of predicted values from both left and
#'   right \code{lambda} indices.
#' @aliases predict.signpen predict.signlog predict.signls
#' @param object fitted \code{\link{signpen}} model object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#'   made. Must be a matrix.
#' @param s value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the entire sequence used to create the model.
#' @param type type of prediction required: \itemize{ \item Type \code{"link"},
#'   for regression it returns the fitted response; for classification it gives
#'   the linear predictors.  \item Type \code{"class"}, only valid for
#'   classification, it produces the predicted class label corresponding to the
#'   maximum probability.}
#' 
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on \code{type}.
#' @keywords perceptron penalty sign consistency
#' @export

predict.signpen = function(object, newx, s = NULL, 
    type = c("class", "link"), ...) NextMethod("predict") 

predict.signlog = function(object, newx, s=NULL, type=c("class", "link"), ...) {
  type = match.arg(type)
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[ , lamlist$left, drop=FALSE] %*% 
      Diagonal(x=lamlist$frac) +
      nbeta[ , lamlist$right, drop=FALSE] %*% 
      Diagonal(x=1-lamlist$frac)
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  nfit = as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  switch(type, link=nfit, class=ifelse(nfit > 0, 1, -1))
} 

predict.signls = function(object, newx, s=NULL, 
    type=c("link"), ...) {
  type = match.arg(type)
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[ , lamlist$left, drop=FALSE] %*% 
      Diagonal(x=lamlist$frac) +
      nbeta[ , lamlist$right, drop=FALSE] %*% 
      Diagonal(x=1-lamlist$frac)
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
} 

#' make predictions from a "cv.signpen" object.
#' 
#' This function makes predictions from a cross-validated \code{signpen} model,
#'   using the stored \code{"signpen.fit"} object, and the optimal value chosen
#'   for \code{lambda}.
#' 
#' This function makes it easier to use the results of cross-validation to make
#'   a prediction.
#' 
#' @param object fitted \code{\link{cv.signpen}} object.
#' @param newx matrix of new values for \code{x} at which predictions are to be
#'   made. Must be a matrix. See documentation for \code{predict.signpen}.
#' @param s value(s) of the penalty parameter \code{lambda} at which
#'   predictions are required. Default is the value \code{s="lambda.1se"} stored
#'   on the CV object. Alternatively \code{s="lambda.min"} can be used. If
#'   \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
#'   used.
#' @param \dots not used. Other arguments to predict.
#' @return The returned object depends on the \dots{} argument which is passed
#'   on to the \code{\link{predict}} method for \code{\link{signpen}} objects.
#' @keywords perceptron penalty sign consistency
#' @export
predict.cv.signpen = function(object, newx, s=c("lambda.1se", "lambda.min"), ...) {
    if (is.numeric(s)) 
        lambda = s else if (is.character(s)) {
        s = match.arg(s)
        lambda = object[[s]]
    } else stop("Invalid form for s")
    predict(object$signpen.fit, newx, s=lambda, ...)
} 


#' print a signpen object
#' 
#' Print the nonzero group counts at each lambda along the signpen path.
#' 
#' Print the information about the nonzero group counts at each lambda step in
#'   the \code{\link{signpen}} object. The result is a two-column matrix with
#'   columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of
#'   the groups that have nonzero within-group coefficients, the \code{Lambda}
#'   column is the the corresponding lambda.
#' 
#' @param x fitted \code{\link{signpen}} object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @return a two-column matrix, the first columns is the number of nonzero
#'   group counts and the second column is \code{Lambda}.
#' @keywords perceptron penalty sign consistency
#' @export
print.signpen = function(x, digits=max(3, getOption("digits") - 3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}


