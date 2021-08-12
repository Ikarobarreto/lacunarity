#' Moments generating Function
#'
#' Calculates the q-moment of frequency distribution of the box masses
#' @param mat is a two-column matrix with x and p(x)
#' @param q is the moment order
#' @return q-moment of a frequency distribution of the box masses

zqs<-function(mat,q){
  Q<-mat$freq/sum(mat$freq)
  z<-(mat$x^q)*Q
  z[is.infinite(z)]<-NA
  zqs<-sum(z,na.rm=TRUE)
  return(zqs)
}
#' lacunarity index
#'
#' Calculates the lacunarity index of a binary series
#' @param x is a vector of 0's and 1's
#' @return y is the lacunarity scaling exponent
#' @return Ds is the lacunarity for the box scale
#' @return s is the box scale
#' @import zoo
#' @import plyr
#' @import stats
#' @examples
#' x<-rbinom(1000,1,0.85)
#' lac(x)
#' @export
#'
lac<-function(x){
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

#' Generalized lacunarity index
#'
#' Calculates the Generalized lacunarity index of a binary series
#' @param x is a vector of 0's and 1's
#' @return s is the box scale
#' @return q is the moment order
#' @return y(q) is the lacunarity scaling exponent gamma(q)
#' @return Dqs is the lacunarity for the box scale and q-moment
#' @import zoo
#' @import plyr
#' @import stats
#' @examples
#' x<-rbinom(1000,1,0.85)
#' genlac(x)
#' @export
genlac<-function(x){
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
