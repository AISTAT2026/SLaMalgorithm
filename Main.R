library(mvtnorm)
Datageneration = function(n, beta, corr, Dist = 'mzNormal', link = 'logit')
{
   beta0 <- beta.0  <- beta
   d <- p <- length(beta.0)-5
   sigmax <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d-1), sigmax)}
   if( Dist == 'nzNormal' ){X  <- rmvnorm(n, rep(1, d-1), sigmax)}
   X =cbind(1, X, rbinom(n, 1, 0.3), rbinom(n, 1, 0.3), rbinom(n, 1, 0.3), runif(n, -2, 2), runif(n, -2, 2))
   
   #X = cbind(1, X)
   p <- d + 1
   if(link == 'logit')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta0)))   
   }else if(link == 'mis')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta0) + (0.3 *X[,2]^2 -0.3 *X[,3]^2) )) 
   }else{
      prob = 1 - exp(-exp(X %*% beta0))
   }
   
   Y = rbinom(n,1,prob)
   list(Y = Y, X = X)
}
#############
A.var <- function(est.score, prob, sub.y, est.m, est.den, class.y.true.v, n, A.thres, type = 'F', type1 = 'IPW')
{
   est.prob = 1 - 1 / (1 + exp(est.score))
   yhat = ifelse(est.prob > A.thres, 1, 0)
  
   if(type == 'F')
   {
     if(type1 == 'IPW')
     {
       u1 = yhat * sub.y; u2 = (class.y.true.v + sub.y)/2; est.u1 = 0; est.u2 = 0
     }else if(type1 == 'SEMI'){ 
       u1 = yhat * sub.y; u2 = (class.y.true.v + sub.y)/2; 
       est.u1 =  yhat * est.prob; est.u2 = (class.y.true.v + est.prob)/2
     }
     est.var = sum((prob - 1) * (( u1 - est.u1 - est.m * (u2 - est.u2) )^2)*prob/(n^2 *  est.den^2) * 4)
   }
   if(type == 'TPR')
   {
     if(type1 == 'IPW')
     {
       u1 = yhat * sub.y; u2 = sub.y; est.u1 = 0; est.u2 = 0
     }else if(type1 == 'SEMI'){ 
       u1 = yhat * sub.y; u2 = sub.y
       est.u1 = yhat * est.prob; est.u2 = est.prob
     }
     est.var = sum((prob - 1) * (( u1 - est.u1 - est.m * (u2 - est.u2) )^2)*prob/(n^2 *  est.den^2))
   }
   if(type == 'PPV')
   {
     if(type1 == 'IPW')
     {
       u1 = yhat * sub.y; u2 = class.y.true.v; est.u1 = 0; est.u2 = 0
       est.var = sum((prob - 1) * (( u1)^2)*prob/(n^2 *  est.den^2))
     }else if(type1 == 'SEMI'){ 
       u1 = yhat * sub.y; u2 = class.y.true.v
       est.u1 = yhat * est.prob; est.u2 = class.y.true.v
       est.var = sum((prob - 1) * (( u1 - est.u1)^2)*prob/(n^2 *  est.den^2))
     }
   } 
   return(est.var)
}
#############
Aopt <- function(full.x, full.y, t.size, sub.x, sub.y, sub.p, batch, c.est, bw, thres, est.M, type ='F', a, b, w = TRUE, setup = FALSE)
{
   if(setup == TRUE)
   {
      pool.x.s = full.x; pool.y.s = full.y;
      sub.y.s = sub.y; sub.x.s = sub.x 
      sub.p.s = sub.p; upd.b.s = c.est; prop.c = NULL; pi.i.s.un=NULL; en.prob1 =NULL
   }else{
      t.n = length(full.y); n = t.size; pool.y.s = full.y;
      est.score = a + b*(full.x %*%  c.est)
      est.prob = 1 - 1 / (1 + exp(est.score))
      yhat = ifelse(est.prob > thres, 1, 0)
      data.eval.en <- data.frame(pool.y=factor(rep(1, length(est.score))),est.score=est.score)
      e.proba.t1 = est.prob
      #e.proba.t1 = 0.9*est.prob + 0.1*0.5
      
      if(type == 'F')
      {
        wegt.e = abs( 2*yhat - 2*(yhat*e.proba.t1) - est.M*(1 - e.proba.t1)) * e.proba.t1 + 
          abs( 2*(yhat*e.proba.t1) - est.M*(e.proba.t1)) *(1 - e.proba.t1)       
        #wegt.e = (yhat - yhat*est.prob) * est.prob + yhat*est.prob *(1 - est.prob)       
        #wegt.e = (est.prob - est.prob*est.prob) * est.prob + est.prob*est.prob *(1 - est.prob)       
        #wegt.e = (1 - est.prob) * est.prob + est.prob *(1 - est.prob)   
        #wegt.e = abs( yhat- est.M) * est.prob 
      }
      if(type == 'TPR')
      {
        wegt.e = abs(yhat - (yhat*e.proba.t1) - est.M*(1 - e.proba.t1)) * e.proba.t1 + 
          abs( (yhat*e.proba.t1) - est.M*(e.proba.t1)) *(1 - e.proba.t1)     
        #wegt.e =  abs(yhat - est.M) * est.prob   
      }
      if(type == 'PPV')
      {
        wegt.e = abs(yhat - (yhat*e.proba.t1)) * e.proba.t1 + 
          abs( (yhat*e.proba.t1)) *(1 - e.proba.t1) 
        #wegt.e = abs(yhat - (est.M*est.prob)) * est.prob + 
        # abs( (est.M*est.prob)) *(1 - est.prob) 
        wegt.e = ifelse(wegt.e == 0, min(wegt.e[wegt.e != 0])/10^10, wegt.e)
      } 
      if(type == 'zero.zero')
      {
        wegt.e = abs( (1-yhat) - ((1-yhat) * est.prob)) * est.prob + 
          abs( ((1-yhat) * est.prob)) *(1 - est.prob) 
        wegt.e = ifelse(wegt.e == 0, min(wegt.e[wegt.e != 0])/10^3, wegt.e)
      } 
      en.prob = wegt.e/sum(wegt.e)

      en.prob1 <- batch * en.prob
      ind.x.en <- (1:n)[runif(n) <= en.prob1]
      
      sub.x.s = rbind(sub.x, full.x[ind.x.en,])
      sub.y.s = c(sub.y, full.y[ind.x.en])
      
      prop.c <- sum(full.y[ind.x.en])
      pool.x.s = full.x[-ind.x.en,]
      pool.y.s = full.y[-ind.x.en]
      pi.i.s.un <- c(sub.p, 1/en.prob1[ind.x.en])
      sub.p.s = c(sub.p, 1/en.prob1[ind.x.en]) + 1
      upd.b.s = c.est
      en.prob1 = c(en.prob1[ind.x.en], en.prob1[-ind.x.en])
   }
   return(list(full.x = pool.x.s, full.y = pool.y.s, sub.x = sub.x.s, pi = pi.i.s.un, en.prob1 = en.prob1,
               sub.y=sub.y.s, sub.p = sub.p.s, c.est = upd.b.s, prop = prop.c, bw = bw, est.M = est.M) )
}
############## entropy####
entropy <- function(full.x, full.y, t.size, sub.x, sub.y, sub.p, batch, c.est, a, b, w = TRUE, setup = FALSE) 
{
   if(setup == TRUE)
   {
      pool.x.en = full.x; pool.y.en = full.y; sub.x.en = sub.x; sub.y.en=sub.y; en.prob1 = NULL
      sub.p.en = sub.p; upd.b.en = c.est; prop.c = NULL; pi.i.en.un = NULL; ind.x.en = NULL
   }else{
      t.n = length(full.y); n = t.size; 
      score = a +  b * full.x %*% c.est
      p_e = c(1 - 1 / (1 + exp(score)))
      wegt.e = p_e* log(p_e) + (1-p_e) * log(1-p_e)
      en.prob = wegt.e/sum(wegt.e)
   
      en.prob1 <- batch * en.prob
      ind.x.en <- (1:n)[runif(n) <= en.prob1]

      sub.x.en = rbind(sub.x, full.x[ind.x.en,])
      sub.y.en = c(sub.y, full.y[ind.x.en])
      
      prop.c <- sum(full.y[ind.x.en])
      pool.x.en = full.x[-ind.x.en,]
      pool.y.en = full.y[-ind.x.en]
      pi.i.en.un <- c(sub.p, 1/en.prob1[ind.x.en])
      sub.p.en = c(sub.p, 1/en.prob1[ind.x.en]) + 1
      upd.b.en = c.est
      en.prob1 = c(en.prob1[ind.x.en], en.prob1[-ind.x.en])
   }   
   return(list(full.x = pool.x.en, full.y = pool.y.en, sub.x = sub.x.en, pi = pi.i.en.un, en.prob1 = en.prob1,
               sub.y=sub.y.en, sub.p = sub.p.en, c.est = upd.b.en, prop = prop.c, ind = ind.x.en) )
}
#############
Exp.Aopt <- function(full.x, full.y, t.size, sub.x, sub.y, sub.p, batch, c.est, bw, thres, est.M, type ='F', a, b, w = TRUE, setup = FALSE)
{
  if(setup == TRUE)
  {
    pool.x.s = full.x; pool.y.s = full.y;
    sub.y.s = sub.y; sub.x.s = sub.x 
    sub.p.s = sub.p; upd.b.s = c.est; prop.c = NULL; pi.i.s.un=NULL; en.prob1 =NULL
  }else{
    t.n = length(full.y); n = t.size; pool.y.s = full.y;
    est.score = a + b*(full.x %*%  c.est)
    est.prob = 1 - 1 / (1 + exp(est.score))
    yhat = ifelse(est.prob > thres, 1, 0)
    data.eval.en <- data.frame(pool.y=factor(rep(1, length(est.score))),est.score=est.score)
    e.proba.t1 = est.prob
    #e.proba.t1 = 0.9*est.prob + 0.1*0.5
    
    if(type == 'F')
    {
      wegt.e = abs( 2*yhat - est.M) * sqrt(e.proba.t1 * (1 - e.proba.t1))
      #wegt.e = abs( 2*yhat - 2*(yhat*e.proba.t1) - est.M*(1 - e.proba.t1)) * e.proba.t1 + 
      #  abs( 2*(yhat*e.proba.t1) - est.M*(e.proba.t1)) *(1 - e.proba.t1)       
      #wegt.e = (yhat - yhat*est.prob) * est.prob + yhat*est.prob *(1 - est.prob)       
      #wegt.e = (est.prob - est.prob*est.prob) * est.prob + est.prob*est.prob *(1 - est.prob)       
      #wegt.e = (1 - est.prob) * est.prob + est.prob *(1 - est.prob)   
      #wegt.e = abs( yhat- est.M) * est.prob 
    }
    if(type == 'PPV')
    {
      wegt.e = abs( (yhat)) *sqrt((1 - e.proba.t1) * e.proba.t1)
      #wegt.e = abs(yhat - (yhat*e.proba.t1)) * e.proba.t1 + 
      #  abs( (yhat*e.proba.t1)) *(1 - e.proba.t1) 
      #wegt.e = abs(yhat - (est.M*est.prob)) * est.prob + 
      # abs( (est.M*est.prob)) *(1 - est.prob) 
      wegt.e = ifelse(wegt.e == 0, min(wegt.e[wegt.e != 0])/10^10, wegt.e)
    } 
    en.prob = wegt.e/sum(wegt.e)
    
    en.prob1 <- batch * en.prob
    ind.x.en <- (1:n)[runif(n) <= en.prob1]
    
    sub.x.s = rbind(sub.x, full.x[ind.x.en,])
    sub.y.s = c(sub.y, full.y[ind.x.en])
    
    prop.c <- sum(full.y[ind.x.en])
    pool.x.s = full.x[-ind.x.en,]
    pool.y.s = full.y[-ind.x.en]
    pi.i.s.un <- c(sub.p, 1/en.prob1[ind.x.en])
    sub.p.s = c(sub.p, 1/en.prob1[ind.x.en]) + 1
    upd.b.s = c.est
    en.prob1 = c(en.prob1[ind.x.en], en.prob1[-ind.x.en])
  }
  return(list(full.x = pool.x.s, full.y = pool.y.s, sub.x = sub.x.s, pi = pi.i.s.un, en.prob1 = en.prob1,
              sub.y=sub.y.s, sub.p = sub.p.s, c.est = upd.b.s, prop = prop.c, bw = bw, est.M = est.M) )
}
###################
Uni <- function(pool.x, pool.y, sub.x, sub.y, t.n, batch, ind) {
   sub.x = rbind(sub.x, pool.x[ind,])
   sub.y = c(sub.y, pool.y[ind])
   pool.x = pool.x[-ind,]
   pool.y = pool.y[-ind]
   upd.b.u <- glm(sub.y~sub.x-1,family = "binomial")$coefficient
   return(list(beta = upd.b.u, Y = y.simp))
}
###################
BIC = function(fit)
{
   tLL <- fit$nulldev - deviance(fit)
   k <- fit$df
   nn <- fit$nobs
   
   BIC<-log(nn)*k - tLL
   #AIC <- -tLL+2*k
   lambda.min = fit$lambda[which.min(BIC)]
   return(lambda.min)
}
#############
A.var.one <- function(est.score, prob, sub.y, n, A.thres, type = 'one.one', type1 = 'IPW')
{
   est.prob = 1 - 1 / (1 + exp(est.score))
   yhat = ifelse(est.prob > A.thres, 1, 0)
   
   if(type == 'one.one')
   {
      if(type1 == 'IPW')
      {
         u1 = yhat * sub.y; est.u1 = 0
      }else if(type1 == 'SEMI'){ 
         u1 = yhat * sub.y; est.u1 =  yhat * est.prob;
      }
      est.var = sum((prob - 1) * ((u1 - est.u1)^2)*prob/(n^2))
   }
   if(type == 'zero.zero')
   {
      if(type1 == 'IPW')
      {
         u1 = (1-yhat) * sub.y; est.u1 = 0
      }else if(type1 == 'SEMI'){ 
         u1 = (1-yhat) * sub.y; est.u1 =  (1-yhat) * est.prob;
      }
      est.var = sum((prob - 1) * ((u1 - est.u1)^2)*prob/(n^2))
   }
   if(type == 'preval')
   {
      if(type1 == 'IPW')
      {
         u1 = sub.y; est.u1 = 0
      }else if(type1 == 'SEMI'){ 
         u1 = sub.y; est.u1 =  est.prob;
      }
      est.var = sum((prob - 1) * ((u1 - est.u1)^2)*prob/(n^2))  
   } 
   return(est.var)
}
##############################################
est.s.prob = function(r, t.score, e.proba.u.res2, band, deg, type = 'Ker')
{
   e.proba.u.res.ker = rep(0,length(e.proba.u.res2));
   if(type == 'Ker')
   {   
      bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (length(r)^(-band))), bandwidth.compute= FALSE)
      data.eval.u <- data.frame(r=factor(rep(1, length(e.proba.u.res2))),t.score=e.proba.u.res2)
      e.proba.u.res.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.u))
   }   
   if(type == 'Spline')
   {   
      reg = glm(r~1+ns(t.score,df=deg),family=binomial)
      e.proba.u.res.ker = 1/predict(reg,newdata=data.frame(t.score = e.proba.u.res2),type="response")
   }
   return(e.proba.u.res.ker)
}   
###############################################
HighData = function(n, beta, corr, Dist = 'mzNormal', link = 'logit')
{
   beta0 <- beta.0  <- beta
   d <- p <- length(beta.0)
   sigmax <- matrix(corr, d-1, d-1) + diag(1-corr, d-1)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d-1), sigmax)}
   X = cbind(1, X)
   p <- d + 1
   if(link == 'logit')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta0)))   
   }else if(link == 'mis')
   {
      prob = 1 - 1 / (1 + exp(c(X %*% beta0) +0.5*X[,3]*X[,4])) 
   }else{
      prob = 1 - exp(-exp(X %*% beta0))
   }
   
   Y = rbinom(n,1,prob)
   list(Y = Y, X = X)
}
###############
getMLE <- function(x, y, w=1) {
   d = ncol(x)
   beta <- rep(0, d)
   dig <- diag(0.00001, d)
   loop  <- 1
   Loop  <- 100
   msg <- "NA"
   while (loop <= Loop) {
      pr <- c(1 - 1 / (1 + exp(x %*% beta)))
      H <- t(x) %*% (pr * (1 - pr) * w * x)
      S <- colSums((y - pr) * w * x)
      tryCatch(
         {shs <- NA
         shs <- solve(H + dig, S) }, 
         error=function(e){
            cat("\n ERROR :", loop, conditionMessage(e), "\n")})
      if (is.na(shs[1])) {
         msg <- "Not converge"
         beta <- loop <- NA
         break
      }
      beta.new <- beta + shs
      tlr  <- sum((beta.new - beta)^2)
      beta  <- beta.new
      if(tlr < 0.000001) {
         msg <- "Successful convergence"
         break
      }
      if (loop == Loop) {
         warning("Maximum iteration reached")
         ## beta <- NA
      }
      loop  <- loop + 1
   }
   list(par=beta, message=msg, iter=loop)
}