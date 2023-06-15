library(pROC);library(MLmetrics);library(PRROC); library(splines);library(np);

source('Main.R')
set.seed(20230401)
tota.b = c(-5.4, rep(0.5,8));re=500; n = 10000; 
A.thres = 0.15; band = 1/2; deg = 10; c1=0.3; c2=0.7
second.est.en.sp <- second.est.en.ker <- second.est.u.sp <- second.est.u.ker <-
  first.est.w.sp <- first.est.w.ker <- a.hat <- b.hat <- c()
auc.t <- t.prauc1 <- propt <- propt.w.en <-propt.w.F <-propt.w.tpr <- propt.w.ppv <-
  propt.u <-  batch.u <- batch.w.en <- batch.w.F <- batch.w.tpr <- batch.w.ppv <- c()
t.measure<- Uni.IPW <- Uni.ker <- Uni.spl <- en.w.IPW <- en.w.ker <- en.w.spl <-
  F.w.IPW <- F.w.ker <- F.w.spl <- tpr.w.IPW <- tpr.w.ker <- tpr.w.spl <-
  ppv.w.IPW <- ppv.w.ker <- ppv.w.spl <- list()

for(i in 1:re)
{
  t.n = 10000
  test.dat = Datageneration(t.n, beta = tota.b,Dist = 'nzNormal', corr = 0.4)
  test.score = c1+c2*(test.dat$X %*% tota.b)
  e.proba.t = 1 - 1 / (1 + exp(test.score))
  test.Y = rbinom(n,1,e.proba.t)
  
  test.pool.y.u = test.Y; test.pool.x.u = test.dat$X
  
  t.inx = sample(1:t.n, t.n/2); 
  test1.pool.x = test.dat$X[t.inx,]; test1.pool.y = test.Y[t.inx]
  test2.pool.x = test.dat$X[-t.inx,]; test2.pool.y = test.Y[-t.inx]
  batch1 = 150; batch2 = 250;
  subt.x.u <- subt.y.u <- subt.p.u <- c()
  set.en <-entropy(test1.pool.x, test1.pool.y, t.n, subt.x.u , subt.y.u, subt.p.u, batch, tota.b, a, b, w = FALSE, setup = TRUE)
  set.en2 <-entropy(test2.pool.x, test2.pool.y, t.n, subt.x.u , subt.y.u, subt.p.u, batch, tota.b, a, b, w = FALSE, setup = TRUE)
  
  #Aopt initial step
  set.F <-Aopt(test2.pool.x, test2.pool.y, t.n, subt.x.u, subt.y.u, subt.p.u, batch, tota.b, bw = NA, A.thres, est.M = NULL, type = 'F', a, b, w = FALSE, setup = TRUE)
  set.tpr <-Aopt(test2.pool.x, test2.pool.y, t.n, subt.x.u, subt.y.u, subt.p.u, batch, tota.b, bw = NA, A.thres, est.M = NULL, type = 'TPR', a, b, w = FALSE, setup = TRUE)
  set.ppv <-Aopt(test2.pool.x, test2.pool.y, t.n, subt.x.u, subt.y.u, subt.p.u, batch, tota.b, bw = NA, A.thres, est.M = NULL, type = 'PPV', a, b, w = FALSE, setup = TRUE)
  
    ##Model Calibration
    set.en = entropy(set.en$full.x, set.en$full.y, t.n/2, set.en$sub.x, set.en$sub.y, set.en$sub.p, batch1, set.en$c.est, 0, 1, w = TRUE)
    length(set.en$sub.y)
    propt.w.en[i] = set.en$prop/length(set.en$sub.y)
    batch.w.en[i] <- length(set.en$sub.y)
    
    est.score.en =  set.en$sub.x %*%  tota.b
    reclib = glm(set.en$sub.y~est.score.en, family = "binomial")$coefficient
    a = reclib[1]; b = reclib[2]; a.hat[i] = reclib[1]; b.hat[i] = reclib[2]
      
    score.en = a + b * est.score.en
    t.proba.re = c(1 - 1 / (1 + exp(score.en)))
    
    e.proba.en.res2 =  (set.en$sub.x %*%  tota.b)
    e.proba.en.res3 =  (set.en$full.x %*%  tota.b) 
    t.score = c(1 - 1 / (1 + exp(set.en$sub.x %*%  tota.b)), 1 - 1 / (1 + exp(set.en$full.x %*%  tota.b)))
    r = c(rep(1,length(e.proba.en.res2)),rep(0,length(e.proba.en.res3)))
    ### Entropy 
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n/2)^(-band)), bandwidth.compute= FALSE)
    data.eval.en <- data.frame(r=factor(rep(1, length(e.proba.en.res2))),t.score=1 - 1 / (1 + exp(set.en$sub.x %*%  tota.b)))
    f.e.proba.en.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.en))

    thres = seq(0,1,0.005)
    yhat.en = t(sapply(t.proba.re, function(x) ifelse(x>thres,1,0)))
    est.proba.t = 1 - 1 / (1 + exp( a+b*(test.dat$X %*% tota.b)))
    class.y.true = t(sapply(est.proba.t, function(x) ifelse(x>thres,1,0)))
    
    deno.w.ker = sum(set.en$sub.y*f.e.proba.en.ker)/sum(f.e.proba.en.ker)
    nume.w.ker = colSums(yhat.en*set.en$sub.y*f.e.proba.en.ker)/sum(f.e.proba.en.ker)
    nume.w.f.ker =  colSums(f.e.proba.en.ker*yhat.en*!set.en$sub.y)/sum(f.e.proba.en.ker)
    FPR.ker = nume.w.f.ker/( 1 - deno.w.ker)
    pre.ind = which(FPR.ker <=A.thres)[1]
    TPR.ker = nume.w.ker[pre.ind]/deno.w.ker
    PPV.ker = nume.w.ker[pre.ind]/colMeans(class.y.true)[pre.ind]; 
    F.ker = 2*TPR.ker* PPV.ker/(TPR.ker + PPV.ker) 

    ##Model evaluation
    a = reclib[1]; b = reclib[2]
    ###Passive
    UniP <- (batch2)/(t.n/2)
    ind.x.u <- (1:(t.n/2))[runif((t.n/2)) <= UniP]
    subt.x.u = test2.pool.x[ind.x.u,]; subt.y.u = test2.pool.y[ind.x.u]
    propt.u[i] <- mean(test2.pool.y[ind.x.u]); subt.p.u = rep(UniP,  length(subt.y.u))
    test.pool.x.u = test2.pool.x[-ind.x.u,]; test.pool.y.u = test2.pool.y[-ind.x.u]
    batch.u[i] <- length(ind.x.u)
    
    ## entropy
    set.en2 = entropy(set.en2$full.x, set.en2$full.y, t.n/2, set.en2$sub.x, set.en2$sub.y, set.en2$sub.p, batch2, set.en2$c.est, a, b, w = TRUE)
    propt.w.en[i] = set.en2$prop/length(set.en2$sub.y)
    batch.w.en[i] <- length(set.en2$sub.y)
    
    ## F
    set.F = Aopt(set.F$full.x, set.F$full.y, length(set.F$full.y), set.F$sub.x, set.F$sub.y, set.F$sub.p, batch2, set.F$c.est, set.F$bw, thres[pre.ind], est.M = F.ker, type = 'F', a, b, w = TRUE)
    propt.w.F[i] = set.F$prop/length(set.F$sub.y)
    batch.w.F[i] <- length(set.F$sub.y)
    
    ## ppv
    set.ppv = Aopt(set.ppv$full.x, set.ppv$full.y, length(set.ppv$full.y), set.ppv$sub.x, set.ppv$sub.y, set.ppv$sub.p, batch2, set.ppv$c.est, set.ppv$bw, thres[pre.ind], est.M = PPV.ker, type = 'PPV', a, b, w = TRUE)
    propt.w.ppv[i] = set.ppv$prop/length(set.ppv$sub.y)
    batch.w.ppv[i] <- length(set.ppv$sub.y)
    
    ####measure
    subt.x.u.t = rbind(set.en$sub.x, subt.x.u); test.pool.x.u.t = rbind(set.en$full.x, test.pool.x.u)
    subt.x.en = rbind(set.en$sub.x, set.en2$sub.x); test.pool.x.en = rbind(set.en$full.x, set.en2$full.x)
    subt.x.F = rbind(set.en$sub.x, set.F$sub.x); test.pool.x.F = rbind(set.en$full.x, set.F$full.x)
    subt.x.ppv = rbind(set.en$sub.x, set.ppv$sub.x); test.pool.x.ppv = rbind(set.en$full.x, set.ppv$full.x)
    u.score = a + b*(subt.x.u.t %*%  tota.b); e.proba.u = 1 - 1 / (1 + exp(u.score))
    u.res.score = a + b*(test.pool.x.u.t %*%  tota.b); e.proba.u.res = 1 - 1 / (1 + exp(u.res.score))
    en.score = a + b*(subt.x.en %*% tota.b); e.proba.en = 1 - 1 / (1 + exp(en.score))
    en.res.score = a + b*(test.pool.x.en %*% tota.b); e.proba.en.res = 1 - 1 / (1 + exp(en.res.score))
    F.score = a + b*(subt.x.F %*% tota.b); e.proba.F = 1 - 1 / (1 + exp(F.score))
    F.res.score = a + b*(test.pool.x.F %*% tota.b); e.proba.F.res = 1 - 1 / (1 + exp(F.res.score))
    ppv.score = a + b*(subt.x.ppv %*% tota.b); e.proba.ppv = 1 - 1 / (1 + exp(ppv.score))
    ppv.res.score = a + b*(test.pool.x.ppv %*% tota.b); e.proba.ppv.res = 1 - 1 / (1 + exp(ppv.res.score))
    
    thres = seq(0,1,0.005)
    est.proba.t = 1 - 1 / (1 + exp( a+b*(test.dat$X %*% tota.b)))
    class.y.true = t(sapply(est.proba.t, function(x) ifelse(x>thres,1,0)))
    class.y.u = t(sapply(e.proba.u, function(x) ifelse(x>thres,1,0)))
    class.y.u.res = t(sapply(e.proba.u.res, function(x) ifelse(x>thres,1,0)))
    class.y.en = t(sapply(e.proba.en, function(x) ifelse(x>thres,1,0)))
    class.y.en.res =  t(sapply(e.proba.en.res, function(x) ifelse(x>thres,1,0)))
    class.y.F = t(sapply(e.proba.F, function(x) ifelse(x>thres,1,0)))
    class.y.F.res =  t(sapply(e.proba.F.res, function(x) ifelse(x>thres,1,0)))
    class.y.ppv = t(sapply(e.proba.ppv, function(x) ifelse(x>thres,1,0)))
    class.y.ppv.res =  t(sapply(e.proba.ppv.res, function(x) ifelse(x>thres,1,0)))
    
    
    ### Uniform-IPW
    subt.p.u = c(set.en$pi, 1/rep(UniP,  length(subt.y.u)))
    subt.y.u.t = c(set.en$sub.y, subt.y.u); 
  
    nume.u = colSums(class.y.u*subt.y.u.t*subt.p.u)/sum(subt.p.u) #colSums(class.y.u*subt.y.u.t*subt.p.u)/(t.n) 
    deno.u = sum(subt.y.u.t*subt.p.u)/sum(subt.p.u) #sum(subt.y.u.t*subt.p.u)/(t.n)
    nume.u.f = colSums(subt.p.u*class.y.u*!subt.y.u.t)/sum(subt.p.u)#colSums(subt.p.u*class.y.u*!subt.y.u.t)/(t.n) 
    
    TPRt.u = nume.u/deno.u
    FPRt.u = nume.u.f/( 1 - deno.u)
    PPVt.u = nume.u/colMeans(class.y.true); PPVt.u = ifelse(PPVt.u>1 , 1,PPVt.u); PPVt.u[colMeans(class.y.true) == 0] = 0
    NPVt.u = (1-FPRt.u) * (1-deno.u)/((1-TPRt.u)*deno.u + (1-FPRt.u)*(1-deno.u))
    Ft.u = 2*TPRt.u* PPVt.u/(TPRt.u + PPVt.u)
    Uni.IPW[[i]] = data.frame(TPRt.u, FPRt.u, PPVt.u, Ft.u, thres, nume.u, deno.u)
    
    ### Uniform-Ker-Estimated prbobability
    e.proba.u.res2 = a+b*(subt.x.u.t %*%  tota.b)
    e.proba.u.res3 = a+b*(test.pool.x.u.t %*%  tota.b)
    t.score = c(1 - 1 / (1 + exp(subt.x.u.t %*%  tota.b)), 1 - 1 / (1 + exp(test.pool.x.u.t %*%  tota.b)))
    r = c(rep(1,length(e.proba.u.res2)),rep(0,length(e.proba.u.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.u <- data.frame(r=factor(rep(1, length(e.proba.u.res2))),t.score=1 - 1 / (1 + exp(subt.x.u.t %*%  tota.b)))
    e.proba.u.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.u))

    e.proba.u.ker = c(e.proba.u.ker)
    nume.u.ker = colSums(class.y.u*subt.y.u.t*e.proba.u.ker)/sum(e.proba.u.ker)#colSums(class.y.u*subt.y.u.t*e.proba.u.ker)/(t.n)
    nume.u.f.ker = colSums(e.proba.u.ker*class.y.u*!subt.y.u.t)/sum(e.proba.u.ker) #colSums(e.proba.u.ker*class.y.u*!subt.y.u.t)/(t.n)
    deno.u.ker = sum(subt.y.u.t*e.proba.u.ker)/sum(e.proba.u.ker)#sum(subt.y.u.t*e.proba.u.ker)/(t.n)
    
    TPR.u.ker = nume.u.ker/deno.u.ker
    FPR.u.ker = nume.u.f.ker/( 1 - deno.u.ker); 
    PPV.u.ker = nume.u.ker/colMeans(class.y.true); PPV.u.ker = ifelse(PPV.u.ker>1 , 1, PPV.u.ker); PPV.u.ker[colMeans(class.y.true) == 0] = 0
    NPV.u.ker = (1-FPR.u.ker)*(sum(e.proba.u.ker*!subt.y.u.t)/t.n)/((1-TPR.u.ker)*deno.u.ker + (1-FPR.u.ker)*(sum(e.proba.u.ker*!subt.y.u.t)/t.n))
    F.u.ker = 2*TPR.u.ker* PPV.u.ker/(TPR.u.ker + PPV.u.ker)
    Uni.ker[[i]] = data.frame(TPR.u.ker, FPR.u.ker, PPV.u.ker, F.u.ker, thres, nume.u.ker, deno.u.ker)
    
    est.score.Avar = c(subt.x.u.t %*% tota.b); class.y.true.v = mean( c(e.proba.u>A.thres, e.proba.u.res>A.thres) )
    AvarSEMIker.u.F[i] = A.var(est.score.Avar, e.proba.u.ker, subt.y.u.t, F.u.ker[thres == A.thres], deno.u.ker + class.y.true.v,class.y.true.v, n, A.thres, type = 'F', type1 = 'SEMI')
    AvarSEMIker.u.ppv[i] = A.var(est.score.Avar, e.proba.u.ker, subt.y.u.t, PPV.u.ker[thres == A.thres], class.y.true.v,class.y.true.v, n, A.thres, type = 'PPV', type1 = 'SEMI')
    
    ### Uniform-Spline-Estimated prbobability
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.u.sp = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(subt.x.u.t %*%  tota.b))),type="response")

    e.proba.u.sp = c(e.proba.u.sp)
    nume.u.ker = colSums(class.y.u*subt.y.u.t*e.proba.u.sp)/sum(e.proba.u.sp) #colSums(class.y.u*subt.y.u.t*e.proba.u.sp)/t.n
    nume.u.f.ker = colSums(e.proba.u.sp*class.y.u*!subt.y.u.t)/sum(e.proba.u.sp) #colSums(e.proba.u.sp*class.y.u*!subt.y.u.t)/t.n
    deno.u.ker = sum(subt.y.u.t*e.proba.u.sp)/sum(e.proba.u.sp) #sum(subt.y.u.t*e.proba.u.sp)/t.n
    
    TPR.u.ker = nume.u.ker/deno.u.ker
    FPR.u.ker = nume.u.f.ker/( 1 - deno.u.ker); 
    PPV.u.ker = nume.u.ker/colMeans(class.y.true); PPV.u.ker = ifelse(PPV.u.ker>1 , 1, PPV.u.ker); PPV.u.ker[colMeans(class.y.true) == 0] = 0
    NPV.u.ker = (1-FPR.u.ker)*(sum(e.proba.u.sp*!subt.y.u.t)/t.n)/((1-TPR.u.ker)*deno.u.ker + (1-FPR.u.ker)*(sum(e.proba.u.sp*!subt.y.u.t)/t.n))
    F.u.ker = 2*TPR.u.ker* PPV.u.ker/(TPR.u.ker + PPV.u.ker)
    Uni.spl[[i]] = data.frame(TPR.u.ker, FPR.u.ker, PPV.u.ker, F.u.ker, thres, nume.u.ker, deno.u.ker)
    
    ### Entropy - IPW
    subt.p.en = c(set.en$pi, set.en2$pi)
    subt.y.en = c(set.en$sub.y, set.en2$sub.y); 
    
    nume.en.w = colSums(class.y.en*subt.y.en*subt.p.en)/sum(subt.p.en) #colSums(class.y.en*subt.y.en*subt.p.en)/t.n
    deno.en.w = sum(subt.y.en*subt.p.en)/sum(subt.p.en) #sum(subt.y.en*subt.p.en)/t.n
    nume.en.f.w = colSums(subt.p.en*class.y.en*!subt.y.en)/sum(subt.p.en) #colSums(subt.p.en*class.y.en*!subt.y.en)/t.n
    TPRt.w.en = nume.en.w/deno.en.w
    FPRt.w.en = nume.en.f.w/( 1 - deno.en.w)
    PPVt.w.en = nume.en.w/colMeans(class.y.true); PPVt.w.en = ifelse(PPVt.w.en>1 , 1, PPVt.w.en); PPVt.w.en[colMeans(class.y.true) == 0] = 0
    NPVt.w.en = (1-FPRt.w.en) * (sum(subt.p.en*!subt.y.en)/t.n)/((1-TPRt.w.en)*deno.en.w + (1-FPRt.w.en)*(sum(subt.p.en*!subt.y.en)/t.n))
    Ft.w.en = 2*TPRt.w.en* PPVt.w.en/(TPRt.w.en + PPVt.w.en)
    en.w.IPW[[i]] = data.frame(TPRt.w.en, FPRt.w.en, PPVt.w.en, Ft.w.en, thres, nume.en.w, deno.en.w)
    
    ### Entropy - Ker
    e.proba.en.res2 = a+b*(subt.x.en %*%  tota.b)
    e.proba.en.res3 = a+b*(test.pool.x.en %*%  tota.b)
    t.score = c(1 - 1 / (1 + exp(subt.x.en %*%  tota.b)), 1 - 1 / (1 + exp(test.pool.x.en %*%  tota.b)))
    r = c(rep(1,length(e.proba.en.res2)),rep(0,length(e.proba.en.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.en <- data.frame(r=factor(rep(1, length(e.proba.en.res2))),t.score=1 - 1 / (1 + exp(subt.x.en %*%  tota.b)))
    e.proba.en.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.en))

    e.proba.en.res.ker = c(e.proba.en.ker)
    deno.w.ker = sum(subt.y.en*e.proba.en.res.ker)/sum(e.proba.en.res.ker) #sum(subt.y.en*e.proba.en.res.ker)/t.n
    nume.w.ker = colSums(class.y.en*subt.y.en*e.proba.en.res.ker)/sum(e.proba.en.res.ker) #colSums(class.y.en*subt.y.en*e.proba.en.res.ker)/t.n
    nume.w.f.ker =  colSums(e.proba.en.res.ker*class.y.en*!subt.y.en)/sum(e.proba.en.res.ker) #colSums(e.proba.en.res.ker*class.y.en*!subt.y.en)/t.n 
    TPR.w.en.ker = nume.w.ker/deno.w.ker
    FPR.w.en.ker = nume.w.f.ker/( 1 - deno.w.ker)
    PPV.w.en.ker = nume.w.ker/colMeans(class.y.true); PPV.w.en.ker = ifelse(PPV.w.en.ker>1 , 1, PPV.w.en.ker); PPV.w.en.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.en.ker = (1-FPR.w.en.ker) * (sum(e.proba.en.res.ker*!subt.y.en)/t.n)/((1-TPR.w.en.ker)*deno.w.ker + (1-FPR.w.en.ker)*(sum(e.proba.en.res.ker*!subt.y.en)/t.n))
    F.w.en.ker = 2*TPR.w.en.ker* PPV.w.en.ker/(TPR.w.en.ker + PPV.w.en.ker)
    en.w.ker[[i]] = data.frame(TPR.w.en.ker, FPR.w.en.ker, PPV.w.en.ker, F.w.en.ker, thres, nume.w.ker, deno.w.ker)
    
    est.score.Avar = c(subt.x.en %*% tota.b); class.y.true.v = mean( c(e.proba.en>A.thres, e.proba.en.res>A.thres) )
    AvarSEMIker.en.F[i] = A.var(est.score.Avar, e.proba.en.res.ker, subt.y.en, F.u.ker[thres == A.thres], deno.w.ker + class.y.true.v,class.y.true.v, n, A.thres, type = 'F', type1 = 'SEMI')
    AvarSEMIker.en.ppv[i] = A.var(est.score.Avar, e.proba.en.res.ker, subt.y.en, PPV.w.en.ker[thres == A.thres], class.y.true.v, class.y.true.v, n, A.thres, type = 'PPV', type1 = 'SEMI')
    
    ### Entropy - Spline
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.en.sp = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(subt.x.en %*%  tota.b))),type="response")
   
    e.proba.en.res.sp = c(e.proba.en.sp)
    deno.w.ker = sum(subt.y.en*e.proba.en.res.sp)/sum(e.proba.en.res.sp) #sum(subt.y.en*e.proba.en.res.sp)/t.n
    nume.w.ker = colSums(class.y.en*subt.y.en*e.proba.en.res.sp)/sum(e.proba.en.res.sp) #colSums(class.y.en*subt.y.en*e.proba.en.res.sp)/t.n
    nume.w.f.ker = colSums(e.proba.en.res.sp*class.y.en*!subt.y.en)/sum(e.proba.en.res.sp) #colSums(e.proba.en.res.sp*class.y.en*!subt.y.en)/t.n 
    TPR.w.en.ker = nume.w.ker/deno.w.ker
    FPR.w.en.ker = nume.w.f.ker/( 1 - deno.w.ker)
    PPV.w.en.ker = nume.w.ker/colMeans(class.y.true); PPV.w.en.ker = ifelse(PPV.w.en.ker>1 , 1, PPV.w.en.ker); PPV.w.en.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.en.ker = (1-FPR.w.en.ker) * (sum(e.proba.en.res.sp*!subt.y.en)/t.n)/((1-TPR.w.en.ker)*deno.w.ker + (1-FPR.w.en.ker)*(sum(e.proba.en.res.sp*!subt.y.en)/t.n))
    F.w.en.ker = 2*TPR.w.en.ker* PPV.w.en.ker/(TPR.w.en.ker + PPV.w.en.ker)
    en.w.spl[[i]] = data.frame(TPR.w.en.ker, FPR.w.en.ker, PPV.w.en.ker, F.w.en.ker, thres, nume.w.ker, deno.w.ker)
    
    ### F - IPW
    subt.p.F = c(set.en$pi, set.F$pi)
    subt.y.F = c(set.en$sub.y, set.F$sub.y); 
    
    nume.F.w = colSums(class.y.F*subt.y.F*subt.p.F)/sum(subt.p.F) #colSums(class.y.F*subt.y.F*subt.p.F)/t.n
    deno.F.w = sum(subt.y.F*subt.p.F)/sum(subt.p.F) #sum(subt.y.F*subt.p.F)/t.n
    nume.F.f.w = colSums(subt.p.F*class.y.F*!subt.y.F)/sum(subt.p.F) #colSums(subt.p.F*class.y.F*!subt.y.F)/t.n 
    TPRt.w.F = nume.F.w/deno.F.w
    FPRt.w.F = nume.F.f.w/( 1 - deno.F.w)
    PPVt.w.F = nume.F.w/colMeans(class.y.true); PPVt.w.F = ifelse(PPVt.w.F>1 , 1, PPVt.w.F); PPVt.w.F[colMeans(class.y.true) == 0] = 0
    NPVt.w.F = (1-FPRt.w.F) * (sum(subt.p.F*!subt.y.F)/t.n)/((1-TPRt.w.F)*deno.F.w + (1-FPRt.w.F)*(sum(subt.p.F*!subt.y.F)/t.n))
    Ft.w.F = 2*TPRt.w.F* PPVt.w.F/(TPRt.w.F + PPVt.w.F)
    F.w.IPW[[i]] = data.frame(TPRt.w.F, FPRt.w.F, PPVt.w.F, Ft.w.F, thres, nume.F.w, deno.F.w)
    
    ### F - Ker
    e.proba.F.res2 = a+b*(subt.x.F%*%  tota.b)
    e.proba.F.res3 = a+b*(test.pool.x.F %*%  tota.b)
    t.score = c(1 - 1 / (1 + exp(subt.x.F %*%  tota.b)), 1 - 1 / (1 + exp(test.pool.x.F %*%  tota.b)))
    r = c(rep(1,length(e.proba.F.res2)),rep(0,length(e.proba.F.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.F <- data.frame(r=factor(rep(1, length(e.proba.F.res2))),t.score=1 - 1 / (1 + exp(subt.x.F %*%  tota.b)))
    e.proba.F.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.F))
 
    e.proba.F.res.ker = c(e.proba.F.ker)
    deno.F.ker = sum(subt.y.F*e.proba.F.res.ker)/sum(e.proba.F.res.ker)
    nume.F.ker = colSums(class.y.F*subt.y.F*e.proba.F.res.ker)/sum(e.proba.F.res.ker)
    nume.F.f.ker = colSums(e.proba.F.res.ker*class.y.F*!subt.y.F) /sum(e.proba.F.res.ker)
    TPR.w.F.ker = nume.F.ker/deno.F.ker
    FPR.w.F.ker = nume.F.f.ker/( 1 - deno.F.ker)
    PPV.w.F.ker = nume.F.ker/colMeans(class.y.true); PPV.w.F.ker = ifelse(PPV.w.F.ker>1 , 1, PPV.w.F.ker); PPV.w.F.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.F.ker = (1-FPR.w.F.ker)*(sum(e.proba.F.res.ker*!subt.y.F)/t.n)/((1-TPR.w.F.ker)*deno.F.ker + (1-FPR.w.F.ker)*(sum(e.proba.F.res.ker*!subt.y.F)/t.n))
    F.w.F.ker = 2*TPR.w.F.ker* PPV.w.F.ker/(TPR.w.F.ker + PPV.w.F.ker)
    F.w.ker[[i]] = data.frame(TPR.w.F.ker, FPR.w.F.ker, PPV.w.F.ker, F.w.F.ker, thres, nume.F.ker, deno.F.ker)
    
    ### F - Spline
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.F.sp = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(subt.x.F %*%  tota.b))),type="response")
    
    e.proba.F.res.sp = c(e.proba.F.sp)
    deno.F.ker = sum(subt.y.F*e.proba.F.res.sp)/sum(e.proba.F.res.sp)
    nume.F.ker = colSums(class.y.F*subt.y.F*e.proba.F.res.sp)/sum(e.proba.F.res.sp)
    nume.F.f.ker = colSums(e.proba.F.res.sp*class.y.F*!subt.y.F) /sum(e.proba.F.res.sp)
    TPR.w.F.ker = nume.F.ker/deno.F.ker
    FPR.w.F.ker = nume.F.f.ker/( 1 - deno.F.ker)
    PPV.w.F.ker = nume.F.ker/colMeans(class.y.true); PPV.w.F.ker = ifelse(PPV.w.F.ker>1 , 1, PPV.w.F.ker); PPV.w.F.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.F.ker = (1-FPR.w.F.ker)*(sum(e.proba.F.res.sp*!subt.y.F)/t.n)/((1-TPR.w.F.ker)*deno.F.ker + (1-FPR.w.F.ker)*(sum(e.proba.F.res.sp*!subt.y.F)/t.n))
    F.w.F.ker = 2*TPR.w.F.ker* PPV.w.F.ker/(TPR.w.F.ker + PPV.w.F.ker)
    F.w.spl[[i]] = data.frame(TPR.w.F.ker, FPR.w.F.ker, PPV.w.F.ker, F.w.F.ker, thres, nume.F.ker, deno.F.ker)

    ### ppv - IPW
    e.proba.ppv.res2 = a+b*(subt.x.ppv %*%  tota.b)
    e.proba.ppv.res3 = a+b*(test.pool.x.ppv %*%  tota.b)
    p.thres = thres[thres <= A.thres]
    e.proba.en.res3 = a+b*(set.en$sub.x %*%  tota.b)
    res3.inx = 1 - 1 / (1 + exp(e.proba.en.res3)) > thres[pre.ind]
    ind.ad = t(sapply(1 - 1 / (1 + exp(e.proba.en.res3)), function(x) ifelse(A.thres>= x & x >p.thres,1,0)))
    
    subt.p.ppv = c(set.en$pi, set.ppv$pi)
    subt.y.ppv = c(set.en$sub.y, set.ppv$sub.y); 
    
    nume.ppv.w = colSums(class.y.ppv*subt.y.ppv*subt.p.ppv)/t.n; nume.ppv.w[thres <= thres[pre.ind]]= colSums(class.y.en[1:length(set.en$sub.y),thres <= thres[pre.ind]]*set.en$sub.y*set.en$pi)/(t.n/2)#nume.ppv.w[thres <= A.thres] = nume.ppv.w[thres <= A.thres] + ad.nume.ppv
    deno.ppv.w = sum(subt.y.ppv*subt.p.ppv)/t.n + sum( (set.en$sub.y*set.en$pi)[!res3.inx] )/t.n#sum(1 - 1 / (1 + exp(e.proba.ppv.res3[!res3.inx])))/(t.n)
    nume.ppv.f.w = colSums(subt.p.ppv*class.y.ppv*!subt.y.ppv)/t.n 
    TPRt.w.ppv = nume.ppv.w/deno.ppv.w; TPRt.w.ppv = ifelse(TPRt.w.ppv>1 , 1, TPRt.w.ppv);
    FPRt.w.ppv = nume.ppv.f.w/( 1 - deno.ppv.w)
    PPVt.w.ppv = nume.ppv.w/colMeans(class.y.true); PPVt.w.ppv = ifelse(PPVt.w.ppv>1 , 1, PPVt.w.ppv); PPVt.w.ppv[colMeans(class.y.true) == 0] = 0
    NPVt.w.ppv = (1-FPRt.w.ppv) * (sum(subt.p.ppv*!subt.y.ppv)/t.n)/((1-TPRt.w.ppv)*deno.ppv.w + (1-FPRt.w.ppv)*(sum(subt.p.ppv*!subt.y.ppv)/t.n))
    Ft.w.ppv = 2*TPRt.w.ppv* PPVt.w.ppv/(TPRt.w.ppv + PPVt.w.ppv)
    ppv.w.IPW[[i]] = data.frame(TPRt.w.ppv, FPRt.w.ppv, PPVt.w.ppv, Ft.w.ppv, thres, nume.ppv.w, deno.ppv.w)
    
    ### ppv - Ker
    t.score = c(1 - 1 / (1 + exp(subt.x.ppv %*%  tota.b)), 1 - 1 / (1 + exp(test.pool.x.ppv %*%  tota.b)))
    r = c(rep(1,length(e.proba.ppv.res2)),rep(0,length(e.proba.ppv.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.ppv <- data.frame(r=factor(rep(1, length(e.proba.ppv.res2))),t.score=1 - 1 / (1 + exp(subt.x.ppv %*%  tota.b)))
    e.proba.ppv.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.ppv))
    
    e.proba.ppv.res.ker = c(e.proba.ppv.ker)
    deno.ppv.ker = sum(subt.y.ppv*e.proba.ppv.res.ker)/t.n # + sum(1 - 1 / (1 + exp(e.proba.ppv.res3[!res3.inx])))/t.n
    nume.ppv.ker = colSums(class.y.ppv*subt.y.ppv*e.proba.ppv.res.ker)/t.n;# nume.ppv.ker[thres <= A.thres] = nume.ppv.ker[thres <= A.thres] + ad.nume.ppv
    nume.ppv.f.ker = colSums(e.proba.ppv.res.ker*class.y.ppv*!subt.y.ppv)/t.n 
    TPR.w.ppv.ker = nume.ppv.ker/deno.ppv.ker; TPR.w.ppv.ker = ifelse(TPR.w.ppv.ker>1 , 1, TPR.w.ppv.ker);
    FPR.w.ppv.ker = nume.ppv.f.ker/( 1 - deno.ppv.ker)
    PPV.w.ppv.ker = nume.ppv.ker/colMeans(class.y.true); PPV.w.ppv.ker = ifelse(PPV.w.ppv.ker>1 , 1, PPV.w.ppv.ker); PPV.w.ppv.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.ppv.ker = (1-FPR.w.ppv.ker)*(sum(e.proba.ppv.res.ker*!subt.y.ppv)/t.n)/((1-TPR.w.ppv.ker)*deno.ppv.ker + (1-FPR.w.ppv.ker)*(sum(e.proba.ppv.res.ker*!subt.y.ppv)/t.n))
    F.w.ppv.ker = 2*TPR.w.ppv.ker* PPV.w.ppv.ker/(TPR.w.ppv.ker + PPV.w.ppv.ker)
    ppv.w.ker[[i]] = data.frame(TPR.w.ppv.ker, FPR.w.ppv.ker, PPV.w.ppv.ker, F.w.ppv.ker, thres, nume.ppv.ker, deno.ppv.ker)
    
    
    ### ppv - Spline
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.ppv.sp = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(subt.x.ppv %*%  tota.b))),type="response")
    
    e.proba.ppv.res.sp = c(e.proba.ppv.sp)
    deno.ppv.ker = sum(subt.y.ppv*e.proba.ppv.res.sp)/t.n# + sum(1 - 1 / (1 + exp(e.proba.ppv.res3[!res3.inx])))/t.n
    nume.ppv.ker = colSums(class.y.ppv*subt.y.ppv*e.proba.ppv.res.sp)/t.n;# nume.ppv.ker[thres <= A.thres] = nume.ppv.ker[thres <= A.thres] + ad.nume.ppv
    nume.ppv.f.ker = colSums(e.proba.ppv.res.sp*class.y.ppv*!subt.y.ppv)/t.n 
    TPR.w.ppv.ker = nume.ppv.ker/deno.ppv.ker; TPR.w.ppv.ker = ifelse(TPR.w.ppv.ker>1 , 1, TPR.w.ppv.ker);
    FPR.w.ppv.ker = nume.ppv.f.ker/( 1 - deno.ppv.ker)
    PPV.w.ppv.ker = nume.ppv.ker/colMeans(class.y.true); PPV.w.ppv.ker = ifelse(PPV.w.ppv.ker>1 , 1, PPV.w.ppv.ker); PPV.w.ppv.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.ppv.ker = (1-FPR.w.ppv.ker)*(sum(e.proba.ppv.res.sp*!subt.y.ppv)/t.n)/((1-TPR.w.ppv.ker)*deno.ppv.ker + (1-FPR.w.ppv.ker)*(sum(e.proba.ppv.res.sp*!subt.y.ppv)/t.n))
    F.w.ppv.ker = 2*TPR.w.ppv.ker* PPV.w.ppv.ker/(TPR.w.ppv.ker + PPV.w.ppv.ker)
    ppv.w.spl[[i]] = data.frame(TPR.w.ppv.ker, FPR.w.ppv.ker, PPV.w.ppv.ker, F.w.ppv.ker, thres, nume.ppv.ker, deno.ppv.ker)
    
    
    print(c(simple_prauc(Uni.IPW[[i]][,1], Uni.IPW[[i]][,2], Uni.IPW[[i]][,3]),
            simple_prauc(Uni.ker[[i]][,1], Uni.ker[[i]][,2] , Uni.ker[[i]][,3]),
            simple_prauc(Uni.spl[[i]][,1], Uni.spl[[i]][,2] , Uni.spl[[i]][,3])))
    print(c(simple_prauc(en.w.IPW[[i]][,1], en.w.IPW[[i]][,2], en.w.IPW[[i]][,3]),
            simple_prauc(en.w.ker[[i]][,1], en.w.ker[[i]][,2] , en.w.ker[[i]][,3]),
            simple_prauc(en.w.spl[[i]][,1], en.w.spl[[i]][,2] , en.w.spl[[i]][,3])))
    print(c(simple_prauc(F.w.IPW[[i]][,1], F.w.IPW[[i]][,2], F.w.IPW[[i]][,3]),
            simple_prauc(F.w.ker[[i]][,1], F.w.ker[[i]][,2] , F.w.ker[[i]][,3]),
            simple_prauc(F.w.spl[[i]][,1], F.w.spl[[i]][,2] , F.w.spl[[i]][,3])))
    print(c(simple_prauc(ppv.w.IPW[[i]][,1], ppv.w.IPW[[i]][,2], ppv.w.IPW[[i]][,3]),
            simple_prauc(ppv.w.ker[[i]][,1], ppv.w.ker[[i]][,2] , ppv.w.ker[[i]][,3]),
            simple_prauc(ppv.w.spl[[i]][,1], ppv.w.spl[[i]][,2] , ppv.w.spl[[i]][,3])))
  
  r.t = suppressMessages(roc(c(test.Y), c(e.proba.t)))
  t.measure[[i]] = as.matrix(coords(r.t, x = "all", input = "threshold", ret = c('tpr', "fpr",'ppv','npv','tp','fp','threshold')));
  auc.t[i] = suppressMessages(auc(test.Y, c(e.proba.t)))
  t.prauc1[i] = PRAUC(e.proba.t, test.Y)
  print(c( t.measure[[i]][which.max(t.measure[[i]][,7]>=A.thres),1] , 
           t.measure[[i]][which.max(t.measure[[i]][,7]>=A.thres),3] ,t.prauc1[i]))
  propt[i] = mean(test.Y)
  print(i)
} 

