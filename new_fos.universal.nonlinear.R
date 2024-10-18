
library(fda)
library(Rcpp)
library(RcppEigen)
sourceCpp("c_fos_universal_nonlinear.cpp")


#########################################
###################################################################
#######################################################################
#' @export
cv.fos.universal.nonlinear=function(X, Y, t.n.basis=20, u.n.basis = 20, K.cv=3, n.u=NULL)
{
  
  n.t=ncol(Y); 
  if(is.null(n.u)){
    n.u=n.t;
  }
  n.sample=nrow(Y) 
  t.y <- seq(0,1,length.out=ncol(Y))
  #shift.x=min(X)
  #scale.x=max(X)-min(X)
  #X=(X-shift.x)/scale.x
  d=scale(X)
  sd_vec=attr(d, "scaled:scale")
  mean_vec=attr(d, "scaled:center")
  X=d
  attr(X, "scaled:center")=NULL
  attr(X, "scaled:scale")=NULL
  attr(X, "dimnames")=NULL
  
  ncomp=ncol(X)+1 
  
  
  bspline.t.obj=create.bspline.basis(c(0,1), t.n.basis, 4)
  bspline.u.obj=create.bspline.basis(c(0,1), u.n.basis, 4)
  
  
  basis.val.t=eval.basis(seq(0, 1, length.out = n.t), bspline.t.obj)
  basis.val.u=eval.basis(seq(0, 1, length.out = n.u), bspline.u.obj)
  
  basis.val.list=list() 
  basis.val.list[[1]]=basis.val.t
  basis.val.list[[2]]=basis.val.u
  
  K.t0=getbasispenalty(bspline.t.obj, 0)
  K.t1=getbasispenalty(bspline.t.obj, 1)
  K.t2=getbasispenalty(bspline.t.obj, 2)
  K.u0=getbasispenalty(bspline.u.obj, 0)
  K.u1=getbasispenalty(bspline.u.obj, 1)
  K.u2=getbasispenalty(bspline.u.obj, 2)
  
  K.u0t2= K.u0 %x% K.t2
  K.u2t0= K.u2 %x% K.t0
  K.u0t0= K.u0 %x% K.t0
  
  x.param.list=list() 
  x.param.list[[1]]=basis.val.t
  x.param.list[[2]]=basis.val.u 
  x.param.list[[3]]=K.t0
  x.param.list[[4]]=K.t2
  x.param.list[[5]]=K.u0
  x.param.list[[6]]=K.u2 
  x.param.list[[7]]=K.u0t2
  x.param.list[[8]]=K.u2t0
  tmp=cbind(  c(1e-1,  1e-1),
              c(1e-2,  1e-2),
              c(1e-2,  1e-3),
              c(1e-3,  1e-2),
              c(1e-3,  1e-3),
              c(1e-3,  1e-5),
              c(1e-5,  1e-5),
              c(1e-5,  1e-7),
              c(1e-6,  1e-4),
              c(1e-6,  1e-6),
              c(1e-6,  1e-8),
              c(1e-7,  1e-9),
              c(1e-7,  1e-7),
              c(1e-7,  1e-5),
              c(1e-8,  1e-6),
              c(1e-8,  1e-8),
              c(1e-8,  1e-10),
              c(1e-9,  1e-10),
              c(1e-9,  1e-9))
  tmp=tmp/n.sample
  x.param.list[[9]]=tmp
  x.param.list[[10]]=mean_vec
  x.param.list[[11]]=sd_vec
  
  X.incr=cbind(rep(1, nrow(X)), X) 
  ind.list=split(sample(n.sample),rep(1:K.cv,length=n.sample))
  X.cv.train=list(K.cv)
  X.cv.valid=list(K.cv)
  Y.cv.train=list(K.cv)
  Y.cv.valid=list(K.cv)
  for(i in 1:K.cv)
  {
    X.cv.train[[i]]=X.incr[-ind.list[[i]], ] 
    X.cv.valid[[i]]=X.incr[ind.list[[i]], ]
    
    Y.cv.train[[i]]=Y[-ind.list[[i]], ] 
    Y.cv.valid[[i]]=Y[ind.list[[i]], ]
  }
  
  fit_cv=cv_universal_nonlinear(X.incr, Y, X.cv.train, X.cv.valid, Y.cv.train, Y.cv.valid, x.param.list)
  
  A=fit_cv$A
  BK=fit_cv$BK
  opt_lambdas=fit_cv$opt_lambdas*n.sample
  
  return(list(A=A, BK=BK, opt.lambdas=opt_lambdas,  fit_cv=fit_cv, x.param.list=x.param.list))
}


######################################
#' @export
pred.fos.universal.nonlinear<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.1<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,1] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.2<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,2] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.3<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,3] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.4<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,4] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.5<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,5] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.6<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,6] 
  return(Yhat)
}


######################################
#' @export
pred.fos.universal.nonlinear.7<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,7] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.8<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,8] 
  return(Yhat)
}


######################################
#' @export
pred.fos.universal.nonlinear.9<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,9] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.10<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,10] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.11<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,11] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.12<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,12] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.13<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,13] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.14<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,14] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.15<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,15] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.16<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,16] 
  return(Yhat)
}


######################################
#' @export
pred.fos.universal.nonlinear.17<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,17] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.18<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,18] 
  return(Yhat)
}


######################################
#' @export
pred.fos.universal.nonlinear.19<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,19] 
  return(Yhat)
}

######################################
#' @export
pred.fos.universal.nonlinear.20<- function(fit.obj,  X.test) #, Y.test)
{
  A=fit.obj$A 
  BK=fit.obj$BK 
  x.param.list=fit.obj$x.param.list 
  basis.val.t=x.param.list[[1]] 
  basis.val.u=x.param.list[[2]] 
  mean_vec=x.param.list[[10]]
  sd_vec=x.param.list[[11]]
  dd=scale(X.test, center=mean_vec, scale=sd_vec)
  X2=cbind(rep(1, nrow(dd)), dd) 
  Y_pred=c_predict(X2,  A,  BK, basis.val.t,  basis.val.u)
  
  Yhat=Y_pred$Yhat[,20] 
  return(Yhat)
}




