#' Epanechnikov Kernel Function
#'
#' Computes the Epanechnikov kernel weights for a given numeric input.
#'
#' The Epanechnikov kernel is defined as
#' \deqn{
#' K(u) = \frac{3}{4}(1-u^2)\mathbf{1}(|u|\le 1)
#' }
#'
#' @param u A numeric vector, matrix, or array of evaluation points.
#'
#' @return A numeric object of the same dimension as \code{u} containing
#' kernel weights.
#'
Ker.epa<-function (u) {
  ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)
}


#' Nadaraya-Watson Smoothing Matrix
#'
#' Constructs a Nadaraya-Watson kernel smoothing matrix based on a distance
#' matrix or index vector.
#'
#' @param tt A vector of indices or a symmetric distance matrix.
#' @param h Bandwidth parameter. If \code{NULL}, an adaptive bandwidth is
#' selected using the empirical 15\% quantile of \code{tt}.
#' @param Ker Kernel function. Default is \code{Ker.epa}.
#' @param w Optional vector of observation weights.
#' @param cv Logical. If \code{TRUE}, diagonal elements are excluded for
#' leave-one-out cross-validation.
#'
#' @return A smoothing matrix with rows normalized to sum to one.
#'
#' @details
#' The function computes kernel weights
#' \deqn{
#' S_{ij} = \frac{K(d_{ij}/h)}{\sum_j K(d_{ij}/h)}
#' }
#' where \eqn{d_{ij}} denotes pairwise distances.
#' @importFrom stats quantile
S.NW<-function (tt, h = NULL, Ker = Ker.epa, w = NULL, cv = FALSE) {
  if (is.matrix(tt)) {
    if (ncol(tt) != nrow(tt)) {
      if (ncol(tt) == 1) {
        tt = as.vector(tt)
        tt = abs(outer(tt, tt, "-"))
      }
    }
  }
  else if (is.vector(tt))
    tt = abs(outer(tt, tt, "-"))
  else stop("Error: incorrect arguments passed")
  if (is.null(h)) {
    h = quantile(tt, probs = 0.15, na.rm = TRUE)
    while (h == 0) {
      h = quantile(tt, probs = pp, na.rm = TRUE)
      pp <- pp + 0.05
    }
  }
  if (cv)
    diag(tt) = Inf
  tt2 <- data.matrix(sweep(tt, 1, h, FUN = "/"))
  k <- Ker(tt2)
  if (is.null(w))
    w <- rep(1, len = ncol(tt))
  k1 <- sweep(k, 2, w, FUN = "*")
  rw <- rowSums(k1, na.rm = TRUE)
  rw[rw == 0] <- 1e-28
  S = k1/rw
  return(S)
}

#' Generalized Cross-Validation Criterion for non-parametric smoothing
#'
#' Computes the generalized cross-validation (GCV) score for a smoothing matrix.
#'
#' @param y A numeric vector or matrix of observations.
#' @param S Nadaraya-Watson Smoothing Matrix.
#'
#' @return A numeric scalar containing the GCV score.
#' The effective degrees of freedom are stored as the attribute \code{"df"}.
#'
#' @details
#' The GCV criterion is defined as
#' \deqn{
#' \mathrm{GCV} =
#' \frac{\|y - Sy\|^2 / n}{(1 - 2\bar{d})}
#' }
#' where \eqn{\bar{d}} is the average diagonal element of \eqn{S}.
#'
GCV.S<-function (y, S){

  n = ncol(S)
  l = 1:n


  if (is.matrix(y) && (ncol(y) > 1)) {
    y2 <- y
  } else if (is.vector(y)) {
    y2 <- y
  } else stop("y is not a vector or matrix")
  # print(dim(S))
  # print(dim(y2))
  y.est = S %*% y2
  # object_size(y2)
  e = y2 - y.est
  res = sum(diag(t(e) %*% e))

  d <- diag(S)[l]
  df <- sum(d)
  if (mean(d, na.rm = TRUE) > 0.5)
    vv = Inf
  else vv = 1/(1 - 2*mean(d, na.rm = TRUE))

  out <- res * vv/n
  attr(out, "df") <- df
  return(out)
}

#' Automatic Selection of Temporal Smoothing Parameter
#'
#' Selects the smoothing bandwidth parameter \code{tau} for functional
#' panel data.
#'
#' @param fts A matrix of functional time series observations with rows
#' corresponding to time points,  columns corresponding to functions, slices corresponding to dimensions.
#' @param p Dimension of the functional observations.
#' @param sn Maximum lag order used in covariance estimation.
#'
#' @return A list containing:
#' \describe{
#'   \item{tau}{Selected smoothing bandwidth.}
#'   \item{res}{Smoothed functional observations.}
#'   \item{n}{Adjusted sample size after trimming.}
#' }
#'
#' @details
#' The bandwidth is selected by minimizing a generalized cross-validation
#' criterion using Brent optimization.
#' @importFrom stats optim quantile
select_tau<-function(fts, p, sn = 6){

  n = nrow(fts)
  N = ncol(fts)
  t = (1:n)/n
  bw = 0.2

  ftsmatrix = apply(fts,1,c)
  ftsgcv <- function(bw){
    S1 <- S.NW(1:ncol(ftsmatrix), floor(bw*n) ,Ker = Ker.epa)
    gcv1 <- GCV.S(t(ftsmatrix), S1)
    return(gcv1[1])
  }
  gcvopt <- optim(0.05, ftsgcv, lower=0.01, upper=max(bw, 0.1), method = "Brent")
  bw = gcvopt$par[1]
  res = locLinSmootherC(bw, t, fts)
  if(p == 1){
    res = array(res[ceiling(n*bw):(n-ceiling(n*bw)),,], dim = c(n-2*ceiling(n*bw)+1, N, p))
  }else{
    res = res[ceiling(n*bw):(n-ceiling(n*bw)),,]
  }

  n = n - 2*ceiling(n*bw) + 1
  ######### select tau ########
  ftsvec = matrix(0, nrow = n, ncol = p^2)
  if(p <= 3){
    ftsvec = matrix(0, nrow = n, ncol = p^2*N*sn)
    for(j in 1:N){
      for(i in (sn+1):n){
        for(k in 1:sn){
          ftsvec[i,(p^2*N*(k-1)+1):(p^2*N*k)] =  as.vector(res[i,j,] %*% t(res[i-k,j,]))
        }
      }
    }
    ftsvec = ftsvec[(sn+1):n,]
  }else{
    for(j in 1:N){
      for(i in (sn+1):n){
        for(k in 1:sn){
          ftsvec[i,] = ftsvec[i,] +  as.vector(res[i,j,] %*% t(res[i-k,j,]))
        }
      }
    }
    ftsvec = ftsvec[(sn+1):n,]/(N*sn)
  }

  ftsgcv2 <- function(tau){
    S1 <- S.NW(1:nrow(ftsvec), floor(tau*n) ,Ker = Ker.epa)
    gcv1 <- GCV.S(ftsvec, S1)
    return(gcv1[1])
  }
  gcvopt2 <- optim(0.35, ftsgcv2, lower=min(n^{-2/5}*4, 0.3), upper=min(n^{-2/5}*6, 0.4), method = "Brent")
  tau = gcvopt2$par[1]

  return(list(tau=tau, res = res, n=n))
}

#' Automatic Block Length Selection
#'
#' Selects the block length parameter used in multiplier bootstrap procedures.
#'
#' @param Qmulti Numeric vector of multiplier statistics.
#' @param lb Lower bound for candidate block lengths.
#' @param ub Upper bound for candidate block lengths.
#'
#' @return An integer corresponding to the selected block length.
#'
#' @details
#' The selected block length minimizes the rolling variance criterion.
#' @importFrom RcppRoll roll_var
select_L<-function(Qmulti, lb, ub){
  if(ub > lb || lb > 1){
    Lmetric = c()
    for(L in lb:ub){
      Lmetric = c(Lmetric, mean(RcppRoll::roll_var(Qmulti, L)))
    }
    Lmetric = RcppRoll::roll_var(Lmetric,3, align="center")
    L = ((lb+1):(ub-1))[which(Lmetric == min(Lmetric))]
  }else{
    print("please give [lb, ub], ub> lb. set L as 7")
    L = 7
  }
  return(L)
}


#' Functional Portmanteau Test
#'
#' Performs a portmanteau-type test for serial dependence in functional
#' time series or panel functional data.
#'
#' @param fts three-dimensional arrays of functional observations, n by N by p, where n is the number of vector functions, N is the number of observations in each function and p is the dimension of the vector functions.
#' @param p Dimension of the vector functions.
#' @param B Number of bootstrap replications.
#' @param sn Maximum lag order.
#' @param lb Lower bound for block length selection is lb \eqn{\times n^{1/5}}$ .
#' @param ub Upper bound for block length selection is ub \eqn{\times n^{1/5}}$.
#' @param tau Smoothing bandwidth for the covariance function.
#' @param trim Trimming parameter.
#' @param auto Logical or numeric indicator. If \code{1}, bandwidth and block
#' length are automatically selected.
#'
#' @return A list containing:
#' \describe{
#'   \item{stat}{Observed test statistic.}
#'   \item{dist}{Bootstrap distribution.}
#'   \item{pvalue}{Bootstrap p-value.}
#' }
#'
#' @details
#' The function implements a bootstrap-based portmanteau test using
#' smoothed covariance estimators.
#'
#' @examples
#' \dontrun{
#' result <- fport(fts, p = 1)
#' result$pvalue
#' }
#'
#' @export
fport<-function(fts, p, B = 1000, sn = 6, lb = 2, ub = 7, tau = 0.2, trim = 1, auto = 1){
  n = nrow(fts)
  N = ncol(fts)
  t = (1:n)/n
  lb  = floor(lb*n^(1/5))
  ub  = floor(ub*n^(1/5))
  sn = floor((log(n))^(2)/sn)

  if(auto == 1){
    result_tau = select_tau(fts, p, sn)
    tau = result_tau$tau
    res = result_tau$res
    n = result_tau$n
  }


  resultQ = ftsQ(tau, trim, t, res, sn)

  if(auto == 1){
    Qmulti = Qmultiplier(tau, trim, t, res, sn)
    L = select_L(Qmulti, lb, ub)
  }else{
    L = lb
  }
  # print(paste("b", bw, "tau",tau, "L", L))
  Qs = Qdist(Qmulti, n, tau, N, B, sn, L)

  pvalue = sum(resultQ/sqrt(sn) < Qs)/B
  return(list(stat = resultQ, dist = Qs, pvalue = pvalue))
}

#' Extended Functional Portmanteau Test
#'
#' Performs an extended portmanteau test for functional panel data with
#' optional universal or exponential covariance structures.
#'
#' @param fts three-dimensional arrays of functional observations, n by N by p, where n is the number of vector functions, N is the number of observations in each function and p is the dimension of the vector functions.
#' @param p Dimension of the vector functions.
#' @param B Number of bootstrap replications.
#' @param sn Maximum lag order.
#' @param lb Lower bound for block length selection is lb \eqn{\times n^{1/5}}$ .
#' @param ub Upper bound for block length selection is ub \eqn{\times n^{1/5}}$
#' @param tau Smoothing bandwidth for the covariance function.
#' @param trim Trimming parameter.
#' @param auto Logical or numeric indicator. If \code{1}, parameters are
#' automatically selected.
#'
#' @return A list containing:
#' \describe{
#'   \item{stat}{Observed test statistic.}
#'   \item{dist}{Bootstrap distribution.}
#'   \item{pvalue}{Bootstrap p-value.}
#' }
#'
#' @details
#' This function extends \code{fport()} to accommodate alternative covariance
#' structures and multivariate functional observations.
#'
#' @examples
#' \dontrun{
#' result <- fport_extend(fts, p = 2)
#' result$pvalue
#' }
#'
#' @export
fport_extend<-function(fts, p, B = 1000, sn = 6, lb = 2, ub = 7, tau = 0.2, trim = 1, auto = 1){
  n = nrow(fts)
  N = ncol(fts)
  t = (1:n)/n
  lb  = floor(lb*n^(1/5))
  ub  = floor(ub*n^(1/5))
  sn = floor((log(n))^(2)/sn)



  if(auto == 1){
    result_tau = select_tau(fts, p, sn)
    tau = result_tau$tau
    res = result_tau$res
    n = result_tau$n
  }


  dotmat = compute_dotmat(res, sn)
  resultQ = ftsQcross_universal(tau, trim, t, res, dotmat, sn)
  if(auto == 1){
    Qmulti = Qmultiplierj_uni(tau, floor(N/2), trim, t, res, dotmat, sn)

    L = select_L(Qmulti, lb, ub)
  }else{
    L = lb
  }
  Qs = Qdist_uni(res, dotmat, t, trim, n ,tau, N, B, sn, L)
  pvalue = sum(resultQ/sqrt(sn) < Qs)/B

  return(list(stat = resultQ, dist = Qs, pvalue = pvalue))
}
