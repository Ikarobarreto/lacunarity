#' lacunarity: lacunarity and generalized lacunarity for binary time series
#'
#' Estimators of the gliding-box lacunarity index and of the generalized
#' (multifractal-like) lacunarity for one-dimensional binary series. See the
#' package vignette, \code{vignette("lacunarity")}, for the underlying theory
#' and worked examples.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{lac}}}{ordinary lacunarity index \eqn{\Lambda(s)} and its
#'     scaling exponent.}
#'   \item{\code{\link{genlac}}}{generalized lacunarity \eqn{\Lambda_q(s)} and
#'     the spectrum of exponents \eqn{\gamma(q)}.}
#' }
#'
#' @keywords internal
"_PACKAGE"

# Internal: validate that `x` is a usable binary series. Stops with a clear
# message on the common user mistakes (non-numeric, NA, non-binary, too short).
.check_binary_series <- function(x, arg = "x") {
  if (missing(x) || is.null(x))
    stop("`", arg, "` must be a numeric vector of 0's and 1's.", call. = FALSE)
  if (is.logical(x)) x <- as.numeric(x)
  if (!is.numeric(x) || !is.null(dim(x)))
    stop("`", arg, "` must be a numeric vector (of 0's and 1's), not a ",
         class(x)[1], ".", call. = FALSE)
  if (anyNA(x))
    stop("`", arg, "` contains missing values (NA); remove them first.",
         call. = FALSE)
  if (!all(x %in% c(0, 1)))
    stop("`", arg, "` must contain only 0's and 1's.", call. = FALSE)
  if (length(x) < 4L)
    stop("`", arg, "` is too short; provide at least 4 observations.",
         call. = FALSE)
  invisible(TRUE)
}

#' Moment generating function of the box masses
#'
#' Internal helper that computes the \eqn{q}-th moment
#' \eqn{Z(q) = \sum_m m^q\, Q(m)} of a box-mass frequency distribution, where
#' \eqn{Q(m)} is the relative frequency of boxes carrying mass \eqn{m}.
#'
#' @param mat a two-column table with the distinct masses \code{x} and their
#'   frequencies \code{freq}.
#' @param q the moment order.
#' @return The \eqn{q}-th moment of the box-mass distribution.
#' @keywords internal
zqs<-function(mat,q){
  Q<-mat$freq/sum(mat$freq)
  z<-(mat$x^q)*Q
  z[is.infinite(z)]<-NA
  zqs<-sum(z,na.rm=TRUE)
  return(zqs)
}
#' Lacunarity index of a binary series
#'
#' Computes the gliding-box lacunarity index \eqn{\Lambda(s)} of a binary time
#' series across dyadic scales, together with its scaling exponent.
#'
#' @details
#' A box of size \eqn{s} is slid one observation at a time along the series and
#' its mass \eqn{m} (the number of ones it covers) is recorded. Writing
#' \eqn{Z(q,s)} for the \eqn{q}-th moment of the resulting box-mass
#' distribution, the lacunarity index is
#' \deqn{\Lambda(s) = \frac{Z(2,s)}{Z(1,s)^2}
#'                  = 1 + \frac{\mathrm{Var}(m)}{\mathrm{mean}(m)^2},}
#' so that \eqn{\Lambda(s) \ge 1}, with equality only for a translationally
#' homogeneous pattern. Larger values indicate gappier, more heterogeneous
#' textures. The scaling exponent \code{y} is the slope of
#' \eqn{\log_2 \Lambda(s)} regressed on \eqn{\log_2 s}. Scales are dyadic,
#' \eqn{s = 2^i}, and capped by the longest run of ones.
#'
#' @param x a binary vector of 0's and 1's.
#' @return A list with components:
#'   \describe{
#'     \item{\code{y}}{the lacunarity scaling exponent \eqn{\gamma}.}
#'     \item{\code{Ds}}{the lacunarity \eqn{\Lambda(s)} at each scale.}
#'     \item{\code{s}}{the dyadic box scales \eqn{s = 2^i}.}
#'   }
#' @references
#' Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of random
#' and deterministic fractal sets. \emph{Physical Review A}, 44(6), 3552-3558.
#'
#' Plotnick, R. E., Gardner, R. H., Hargrove, W. W., Prestegaard, K. and
#' Perlmutter, M. (1996). Lacunarity analysis: a general technique for the
#' analysis of spatial patterns. \emph{Physical Review E}, 53(5), 5461-5468.
#' @seealso \code{\link{genlac}} for the generalized lacunarity spectrum.
#' @import zoo
#' @import plyr
#' @import stats
#' @examples
#' x <- rbinom(1000, 1, 0.85)
#' lac(x)
#' @export
#'
lac<-function(x){
  .check_binary_series(x)
  requireNamespace("zoo")
  if(sum(x)==0||sum(x)==length(x)) {
    warning("No gaps where found")
  }
  else {
    mod<-rle(x)
    mod1<-cbind.data.frame(as.vector(mod$lengths),as.vector(mod$values))
    mod2<-mod1[which(mod1[,2]==1),]
    p<-ceiling(log2(max(mod2)))
    ml<-matrix(NA,ncol=p,nrow=length(x))
    r<-matrix(NA,nrow=p,ncol=1)
    for(i in 1:p){
      max<-length(x)-(2^i)+1
      ml[(1:max),i]<-t(rollapply(x,width=(2^i),FUN=sum,by=1))
      r[i]<-2^i
    }
    LAM1<-vector("numeric",length = p)
    for(i in 1:p){
      LAM1[i]<-(zqs(na.omit(count(ml[,i])),q=2)/
                  zqs(na.omit(count(ml[,i])),1)^2)
    }
    beta<-(lm(log2(LAM1)~log2(r)-1)$coeff[[1]])
    lac<-list(beta,LAM1,r)
    names(lac)<-c("y","Ds","s")
    return(lac)
  }
}

#' Generalized lacunarity of a binary series
#'
#' Computes the generalized (multifractal-like) lacunarity \eqn{\Lambda_q(s)} of
#' a binary time series and the spectrum of scaling exponents \eqn{\gamma(q)}.
#'
#' @details
#' The ordinary lacunarity is extended to an arbitrary moment order \eqn{q} by
#' \deqn{\Lambda_q(s) = \left[ \frac{Z(2q,s)}{Z(q,s)^2} \right]^{1/q},}
#' where \eqn{Z(q,s)} is the \eqn{q}-th moment of the gliding-box mass
#' distribution at scale \eqn{s}. Large positive \eqn{q} emphasises dense boxes
#' and negative \eqn{q} emphasises sparse boxes, so the curve \eqn{q \mapsto
#' \gamma(q)}, with \eqn{\gamma(q)} the slope of \eqn{\log_{10}\Lambda_q(s)} on
#' \eqn{\log_{10} s}, describes how gaps of different magnitudes scale. Orders
#' \eqn{q} range over \eqn{\{-10, \dots, 10\} \setminus \{0\}}.
#'
#' @param x a binary vector of 0's and 1's.
#' @return A list with components:
#'   \describe{
#'     \item{\code{s}}{the dyadic box scales \eqn{s = 2^i}.}
#'     \item{\code{q}}{the moment orders.}
#'     \item{\code{yq}}{the generalized scaling exponents \eqn{\gamma(q)}.}
#'     \item{\code{Dqs}}{the matrix of generalized lacunarities
#'       \eqn{\Lambda_q(s)} (rows index \code{q}, columns index \code{s}).}
#'   }
#' @references
#' Vernon-Carter, J., Lobato-Calleros, C., Escarela-Perez, R., Rodriguez, E. and
#' Alvarez-Ramirez, J. (2009). A suggested generalization for the lacunarity
#' index. \emph{Physica A}, 388(20), 4305-4314.
#'
#' Allain, C. and Cloitre, M. (1991). Characterizing the lacunarity of random
#' and deterministic fractal sets. \emph{Physical Review A}, 44(6), 3552-3558.
#' @seealso \code{\link{lac}} for the ordinary lacunarity index.
#' @import zoo
#' @import plyr
#' @import stats
#' @examples
#' x <- rbinom(1000, 1, 0.85)
#' genlac(x)
#' @export
genlac<-function(x){
  .check_binary_series(x)
  requireNamespace("zoo")
  if(sum(x)==0||sum(x)==length(x)) {
    warning("No gaps where found")
  }
  else {
    mod<-rle(x)
    mod1<-cbind.data.frame(as.vector(mod$lengths),as.vector(mod$values))
    mod2<-mod1[which(mod1[,2]==1),]
    p<-ceiling(log2(max(mod2)))
    M1<-1
    sb<-2^(1:p)
    ml<-matrix(NA,ncol=p,nrow=length(x))
    for(i in 1:p){
      max<-length(x)-sb[i]+1
      ml[(1:max),i]<-t(rollsum(x,k=sb[i]))
    }
    Dqs<-matrix(NA,ncol=p,nrow=20)
    l<-seq(from=-10,to=10,length.out = 21)
    l<-l[which(l!=0)]
    for(i in 1:p){
      for(j in 1:20){
        Dqs[j,i]<-(zqs(na.omit(count(ml[,i])),q=2*l[j])/
                     zqs(na.omit(count(ml[,i])),q=l[j])^2)^(1/l[j])
      }
    }
    yk<-vector(,length = 20)
    for(i in 1:20){
      yk[i]<-sum(log10(t(Dqs[i,]))*log10(sb))/sum(log10(t(sb))^2)
    }
    mlac<-list(sb,l,yk,Dqs)
    names(mlac)<-c("s","q","yq","Dqs")
    return(mlac)
  }
}
