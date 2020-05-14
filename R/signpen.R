#' perceptron penalty for sign consistency
#'
#' @param x A matrix with \eqn{N} rows and \eqn{p} columns for predictors.
#' @param y A vector of length \eqn{p} for responses.
#' @param method A character string specifying the loss function to use. 
#'   Options are 
#'     "ls": least squared regression;
#'     "logit": logistic regression.
#'   Default is "ls".
#' @param nlambda The number of \code{lambda} values, 
#'   i.e., length of the \code{lambda} sequence. Default is 100. 
#' @param lambda.factor Takes the value of 0.0001 if \eqn{N >= p} or 0.01 if \eqn{N < p}.
#'   Takes no effect when user specifies a \code{lambda} sequence.
#' @param lambda An optional user-supplied \code{lambda} sequence. 
#'   If \code{lambda = NULL} (default), the program computes 
#'   its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.factor}; 
#'   otherwise, the program uses the user-specified one. 
#'   Since the program will automatically sort user-defined \code{lambda} sequence in decreasing order, 
#'   it is better to supply a decreasing sequence.
#' @param lambda2 The L2 tuning parameter \eqn{\lambda_2}{lambda2}.
#' @param lambda3 The L3 tuning parameter \eqn{\lambda_3}{lambda3}.
#' @param gam The gamma vector of length \eqn{p} for prior information on the signs of \code{beta}.
#' @param pf A vector of length \eqn{p}{p} representing the L1 penalty weights 
#'  \code{pf} can be 0 for some predictor(s), leading to including the predictor(s) all the time. 
#'   One suggested choice of \code{pf} is \eqn{{(\beta + 1/n)}^{-1}}{1/(beta+1/n)}, 
#'   where \eqn{n} is the sample size and \eqn{\beta}{beta} is the coefficents obtained. 
#'   Default is 1 for all predictors (and infinity if some predictors are listed in \code{exclude}).
#' @param pf2 A vector of length \eqn{p}{p} for L2 penalty factor for adaptive L1 or adaptive elastic net. 
#'   To allow different L2 shrinkage, user can set \code{pf2} to be different L2 penalty weights 
#'   for each coefficient of \eqn{\beta}{beta}. \code{pf2} can be 0 for some variables, 
#'   indicating no L2 shrinkage. Default is 1 for all predictors.
#' @param exclude Whether to exclude some predictors from the model. 
#'   This is equivalent to adopting an infinite penalty factor when excluding some predictor. 
#'   Default is none. 
#' @param dfmax Restricts at most how many predictors can be incorporated in the model. Default is \eqn{p+1}. 
#'   This restriction is helpful when \eqn{p} is large, provided that a partial path is acceptable. 
#' @param pmax Restricts the maximum number of variables ever to be nonzero; 
#'   e.g, once some \eqn{\beta} enters the model, it counts once. 
#'   The count will not change when the \eqn{\beta} exits or re-enters the model. 
#'   Default is \code{min(dfmax*1.2,p)}.
#' @param standardize Whether to standardize the data. If \code{TRUE}, 
#'   \code{\link{signpen}} normalizes the predictors such that each column has 
#'   sum squares\eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N} of one. 
#'   Note that x is always centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}) 
#'   no matter \code{standardize} is \code{TRUE} or \code{FALSE}. 
#'   \code{\link{signpen}} always returns coefficient \code{beta} on the original scale.  
#'   Default value is \code{TRUE}.
#' @param eps The algorithm stops when 
#'   (i.e. \eqn{\max_j(\beta_j^{new}-\beta_j^{old})^2}{max(j)(beta_new[j]-beta_old[j])^2} 
#'   is less than \code{eps}, where \eqn{j=0,\ldots, p}. Defaults value is \code{1e-8}.
#' @param maxit Restricts how many outer-loop iterations are allowed. 
#'   Default is 1e6. Consider increasing \code{maxit} when the algorithm does not converge.
#'
#' @details
#' The \code{\link{signpen}} minimizes the sparse penalized loss function with the perceptron penalty, 
#'   \deqn{L(y, X, \beta)/N + \lambda_1||\beta||_1 + 0.5\lambda_2||\beta||_2^2 + \lambda_3 (-\beta_j \gamma_j)_+. }
#'   {L(y, X, beta))/N + lambda1 * ||beta||_1^1 + 0.5 * lambda2 * ||beta||_2^2 + lambda3 (-beta_j gamma_j)_+. }
#'   The values of \code{lambda2} \code{lambda3} are user-specified.
#'   
#'   When the algorithm do not converge or run slow, consider increasing \code{eps}, decreasing
#'   \code{nlambda}, or increasing \code{lambda.factor} before increasing
#'   \code{maxit}.
#' @return
#' An object with S3 class \code{\link{signpen}}.
#'  \item{b0}{A vector of length \code{length(lambda)} representing the intercept at each \code{lambda} value.}
#'  \item{beta}{A matrix of dimension \code{p*length(lambda)} representing the coefficients at each \code{lambda} value. 
#'   The matrix is stored as a sparse matrix  (\code{Matrix} package). 
#'     To convert it into normal type matrix use \code{as.matrix()}.}
#'  \item{df}{The number of nonzero coefficients at each \code{lambda}.}
#'  \item{dim}{The dimension of coefficient matrix, i.e., \code{p*length(lambda)}.}
#'  \item{lambda}{The \code{lambda} sequence that was actually used.}
#'  \item{npasses}{Total number of iterations for all lambda values. }
#'  \item{jerr}{Warnings and errors; 0 if no error.}
#'  \item{call}{The call that produced this object.}
#' 
#' @keywords perceptron penalty sign consistency
#' @useDynLib signpen
#' @export
#'
#-----------------------------------------------------------------------------
signpen = function(x, y, method = c("ls", "logit"), nlambda=100, 
  lambda.factor=ifelse(nobs < nvars, 0.01, 1e-04), 
  lambda=NULL, lambda2=0, lambda3=1, gam=rep(1, nvars),
  pf=rep(1, nvars), pf2=rep(1, nvars), exclude, dfmax=nvars + 1, 
  pmax=min(dfmax * 1.2, nvars), standardize=TRUE, 
  eps=1e-08, maxit=1e+06) {
  ####################################################################
  #data setup
  method = match.arg(method)
  this.call = match.call()
  y = drop(y)
  x = as.matrix(x)
  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  vnames = colnames(x)
  if (is.null(vnames)) 
    vnames = paste0("V", seq(nvars))
  if (length(y) != nobs) 
    stop("x and y have different number of observations")
  ####################################################################
  #parameter setup
  if (length(pf) != nvars) 
    stop("The size of L1 penalty factor must be same as the number of input variables.")
  if (length(pf2) != nvars) 
    stop("The size of L2 penalty factor must be same as the number of input variables.")
  if (lambda2 < 0) 
    stop("lambda2 must be non-negative")
  if (lambda3 < 0) 
    stop("lambda2 must be non-negative")
  if (length(gam) != nvars) 
    stop("The size of gam vector must be same as the number of input variables.")
  maxit = as.integer(maxit)
  lam2 = as.double(lambda2)
  lam3 = as.double(lambda3)
  gam = as.double(gam)
  pf = as.double(pf)
  pf2 = as.double(pf2)
  isd = as.integer(standardize)
  eps = as.double(eps)
  dfmax = as.integer(dfmax)
  pmax = as.integer(pmax)
  if (!missing(exclude)) {
    jd = match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) 
      stop("Some excluded variables are out of range.")
    jd = as.integer(c(length(jd), jd))
  } else jd = as.integer(0)
  ####################################################################
  #lambda setup
  nlam = as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1.")
    flmin = as.double(lambda.factor)
    ulam = double(1)
  } else {
    #flmin=1 if user define lambda
    flmin = as.double(1)
    if (any(lambda < 0)) 
      stop("The values of lambda should be non-negative.")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  ####################################################################
  fit = switch(method, 
    ls = signls(x, y, nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, 
      pf, pf2, maxit, lam2, lam3, gam, nobs, nvars, vnames), 
    logit = signlog(x, y, nlam, flmin, ulam, isd, eps, dfmax, pmax, jd, 
      pf, pf2, maxit, lam2, lam3, gam, nobs, nvars, vnames), 
    )
  if (is.null(lambda)) 
    fit$lambda = lamfix(fit$lambda)
  fit$call = this.call
  ####################################################################
  class(fit) = c("signpen", class(fit))
  fit
} 

signls = function(x, y, nlam, flmin, ulam, isd, eps, dfmax, pmax, 
  jd, pf, pf2, maxit, lam2, lam3, gam, nobs, nvars, vnames) {
  ####################################################################
  #data setup
  y = as.numeric(y)
  ####################################################################
  # call Fortran core
  fit = .Fortran("ls_sign", lam2, lam3, gam, nobs, nvars, 
    as.double(x), as.double(y), jd, pf, pf2, dfmax, 
    pmax, nlam, flmin, ulam, eps, isd, maxit, 
    nalam=integer(1), b0=double(nlam), beta=double(pmax * nlam), 
    ibeta=integer(pmax), nbeta=integer(nlam), alam=double(nlam), 
    npass=integer(1), jerr=integer(1))
  #################################################################
  # output
  outlist = getoutput(fit, maxit, pmax, nvars, vnames)
  outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) = c("signls")
  outlist
} 


signlog = function(x, y, nlam, flmin, ulam, isd, eps, dfmax, pmax, 
  jd, pf, pf2, maxit, lam2, lam3, gam, nobs, nvars, vnames) {
  ####################################################################
  #data setup
  y = as.factor(y)
  y = c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  ####################################################################
  # call Fortran core
  fit = .Fortran("log_sign", lam2, lam3, gam, nobs, nvars, 
    as.double(x), as.double(y), jd, pf, pf2, dfmax, 
    pmax, nlam, flmin, ulam, eps, isd, maxit, 
    nalam=integer(1), b0=double(nlam), beta=double(pmax * nlam), 
    ibeta=integer(pmax), nbeta=integer(nlam), alam=double(nlam), 
    npass=integer(1), jerr=integer(1))
  #################################################################
  # output
  outlist = getoutput(fit, maxit, pmax, nvars, vnames)
  outlist = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) = c("signlog")
  outlist
} 
