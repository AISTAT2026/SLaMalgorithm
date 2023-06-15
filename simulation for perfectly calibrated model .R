library(pROC); library(PRROC); library(splines);library(np);library(MLmetrics);

source('Main.R')
set.seed(20230401)
tota.b = c(-4.3, rep(0.5,8));re=500; n = 10000; 
 A.thres = 0.15; band = 1/3; deg = 8; a=0; b=1
auc.t <- t.prauc1 <- propt <- propt.w.en <-propt.w.F <-propt.w.tpr <- propt.w.ppv <-
  propt.u <-  batch.u <- batch.w.en <- batch.w.F <- batch.w.tpr <- batch.w.ppv <- c()
t.measure<- Uni.IPW <- Uni.ker <- Uni.spl <- en.w.IPW <- en.w.ker <- en.w.spl <-
  F.w.IPW <- F.w.ker <- F.w.spl <- tpr.w.IPW <- tpr.w.ker <- tpr.w.spl <-
  ppv.w.IPW <- ppv.w.ker <- ppv.w.spl <- list()
AvarIPW.u.ppv <- AvarIPW.en.ppv<- AvarIPW.u.F <- AvarIPW.en.F<-
  AvarSEMIker.u.ppv <- AvarSEMIker.en.ppv <-  AvarSEMIker.u.F <- AvarSEMIker.en.F <-
   AvarSEMIspl.u.ppv <- AvarSEMIspl.en.ppv <-AvarSEMIspl.u.F <- AvarSEMIspl.en.F <- c()

for(i in 1:re)
{
  test.dat = Datageneration(n, beta = tota.b, Dist = 'nzNormal', corr = 0.4)
  test.pool.x = test.dat$X; test.pool.y = test.dat$Y
  batch = 350;
  n1 = length(test.pool.y)
  test.score = a+b*(test.pool.x %*% tota.b)
  e.proba.t = 1 - 1 / (1 + exp(test.score))
  test.pool.y.u = test.pool.y; test.pool.x.u = test.pool.x
  subt.x.u <- subt.y.u <- subt.p.u <- c()
  set.en <-entropy(test.pool.x, test.pool.y, n1, subt.x.u , subt.y.u, subt.p.u, batch, tota.b, a, b, w = FALSE, setup = TRUE)
  
    t.n = length(test.pool.y.u) 
    class.y.true.v = colMeans(e.proba.t>A.thres)
    
    ###Passive
    UniP <- (batch)/t.n
    ind.x.u <- (1:t.n)[runif(t.n) <= UniP]
    subt.x.u = rbind(subt.x.u, test.pool.x.u[ind.x.u,]);subt.y.u = c(subt.y.u, test.pool.y.u[ind.x.u])
    propt.u[i] <- mean(test.pool.y[ind.x.u]); subt.p.u = c(subt.p.u, rep(UniP,  length(subt.y.u)))
    test.pool.x.u = test.pool.x.u[-ind.x.u,]; test.pool.y.u = test.pool.y.u[-ind.x.u]
    batch.u[i] <- length(subt.y.u)
    
    ## entropy
    set.en = entropy(set.en$full.x, set.en$full.y, n1, set.en$sub.x, set.en$sub.y, set.en$sub.p, batch, set.en$c.est, a, b, w = TRUE)
    propt.w.en[i] = set.en$prop/length(set.en$sub.y)
    batch.w.en[i] <- length(set.en$sub.y)

    ####measure
    u.score = a + b*(subt.x.u %*%  tota.b); e.proba.u = 1 - 1 / (1 + exp(u.score))
    u.res.score = a + b*(test.pool.x.u %*%  tota.b); e.proba.u.res = 1 - 1 / (1 + exp(u.res.score))
    en.score = a + b*(set.en$sub.x %*% tota.b); e.proba.en = 1 - 1 / (1 + exp(en.score))
    en.res.score = a + b*(set.en$full.x %*% tota.b); e.proba.en.res = 1 - 1 / (1 + exp(en.res.score))
 
    thres = seq(0,1,0.0005)
    class.y.true = t(sapply(e.proba.t, function(x) ifelse(x>thres,1,0)))
    class.y.u = t(sapply(e.proba.u, function(x) ifelse(x>thres,1,0)))
    class.y.u.res = t(sapply(e.proba.u.res, function(x) ifelse(x>thres,1,0)))
    class.y.en = t(sapply(e.proba.en, function(x) ifelse(x>thres,1,0)))
    class.y.en.res =  t(sapply(e.proba.en.res, function(x) ifelse(x>thres,1,0)))

    
    ### Uniform-IPW
    nume.u = colMeans(class.y.u*subt.y.u); deno.u = mean(subt.y.u) 
    nume.u.f = colMeans(class.y.u*!subt.y.u) #colMeans(class.y.true)-nume.u
    
    TPRt.u = nume.u/deno.u
    FPRt.u = nume.u.f/( 1 - deno.u)
    PPVt.u = nume.u/colMeans(class.y.true); PPVt.u = ifelse(PPVt.u>1 , 1,PPVt.u); PPVt.u[colMeans(class.y.true) == 0] = 0
    NPVt.u = (1-FPRt.u) * (1-deno.u)/((1-TPRt.u)*deno.u + (1-FPRt.u)*(1-deno.u))
    Ft.u = 2*TPRt.u* PPVt.u/(TPRt.u + PPVt.u)
    Uni.IPW[[i]] = data.frame(TPRt.u, FPRt.u, PPVt.u, Ft.u, thres, nume.u, deno.u)
    
    ### Uniform-Ker-Estimated prbobability
    e.proba.u.res2 = a+b*(subt.x.u %*%  tota.b)
    e.proba.u.res3 = a+b*(test.pool.x.u %*%  tota.b)
    t.score = c(1 - 1 / (1 + exp(subt.x.u %*%  tota.b)), 1 - 1 / (1 + exp(test.pool.x.u %*%  tota.b)))
    r = c(rep(1,length(e.proba.u.res2)),rep(0,length(e.proba.u.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.u <- data.frame(r=factor(rep(1, length(e.proba.u.res2))),t.score=1 - 1 / (1 + exp(subt.x.u %*%  tota.b)))
    e.proba.u.res.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.u))
    
    nume.u.ker = colSums(class.y.u*subt.y.u*e.proba.u.res.ker)/n1
    nume.u.f.ker = colSums(e.proba.u.res.ker*class.y.u*!subt.y.u)/n1  #colMeans(class.y.true)-nume.u
    deno.u.ker = sum(subt.y.u*e.proba.u.res.ker)/n1
    
    TPR.u.ker = nume.u.ker/deno.u.ker
    FPR.u.ker = nume.u.f.ker/( 1 - deno.u.ker); 
    PPV.u.ker = nume.u.ker/colMeans(class.y.true); PPV.u.ker = ifelse(PPV.u.ker>1 , 1, PPV.u.ker); PPV.u.ker[colMeans(class.y.true) == 0] = 0
    NPV.u.ker = (1-FPR.u.ker)*(sum(e.proba.u.res.ker*!subt.y.u)/n1)/((1-TPR.u.ker)*deno.u.ker + (1-FPR.u.ker)*(sum(e.proba.u.res.ker*!subt.y.u)/n1))
    F.w.ker = 2*TPR.u.ker* PPV.u.ker/(TPR.u.ker + PPV.u.ker)
    Uni.ker[[i]] = data.frame(TPR.u.ker, FPR.u.ker, PPV.u.ker, F.w.ker, thres, nume.u.ker, deno.u.ker)
    
    est.score.Avar = c(subt.x.u %*% tota.b)
    pre.ind = which(FPR.u.ker <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarSEMIker.u.F[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, F.w.ker[pre.ind], deno.u.ker + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'SEMI')
    AvarSEMIker.u.ppv[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, PPV.u.ker[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'SEMI')
    
    ### Uniform-Spline-Estimated prbobability
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.u.res.ker = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(subt.x.u %*%  tota.b))),type="response")
    
    nume.u.ker = colSums(class.y.u*subt.y.u*e.proba.u.res.ker)/n1
    nume.u.f.ker = colSums(e.proba.u.res.ker*class.y.u*!subt.y.u)/n1 
    deno.u.ker = sum(subt.y.u*e.proba.u.res.ker)/n1
    
    TPR.u.ker = nume.u.ker/deno.u.ker
    FPR.u.ker = nume.u.f.ker/( 1 - deno.u.ker); 
    PPV.u.ker = nume.u.ker/colMeans(class.y.true); PPV.u.ker = ifelse(PPV.u.ker>1 , 1, PPV.u.ker); PPV.u.ker[colMeans(class.y.true) == 0] = 0
    NPV.u.ker = (1-FPR.u.ker)*(sum(e.proba.u.res.ker*!subt.y.u)/n1)/((1-TPR.u.ker)*deno.u.ker + (1-FPR.u.ker)*(sum(e.proba.u.res.ker*!subt.y.u)/n1))
    F.w.ker = 2*TPR.u.ker* PPV.u.ker/(TPR.u.ker + PPV.u.ker)
    Uni.spl[[i]] = data.frame(TPR.u.ker, FPR.u.ker, PPV.u.ker, F.w.ker, thres, nume.u.ker, deno.u.ker)
    
    est.score.Avar = c(subt.x.u %*% tota.b)
    pre.ind = which(FPR.u.ker <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarSEMIspl.u.F[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, F.w.ker[pre.ind], deno.u.ker + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'SEMI')
    AvarSEMIspl.u.ppv[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, PPV.u.ker[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'SEMI')
    
    est.score.Avar = c(subt.x.u %*% tota.b)
    pre.ind = which(FPRt.u <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarIPW.u.F[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, Ft.u[pre.ind], deno.u + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'IPW')
    AvarIPW.u.ppv[i] = A.var(est.score.Avar, 1/subt.p.u, subt.y.u, PPVt.u[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'IPW')
    
    ### Entropy - IPW
    nume.en.w = colSums(class.y.en*set.en$sub.y*set.en$pi)/n1
    deno.en.w = sum(set.en$sub.y*set.en$pi)/n1
    nume.en.f.w = colSums(set.en$pi*class.y.en*!set.en$sub.y)/n1
    TPRt.w.en = nume.en.w/deno.en.w
    FPRt.w.en = nume.en.f.w/( 1 - deno.en.w)
    PPVt.w.en = nume.en.w/colMeans(class.y.true); PPVt.w.en = ifelse(PPVt.w.en>1 , 1, PPVt.w.en); PPVt.w.en[colMeans(class.y.true) == 0] = 0
    NPVt.w.en = (1-FPRt.w.en) * (sum(set.en$pi*!set.en$sub.y)/n1)/((1-TPRt.w.en)*deno.en.w + (1-FPRt.w.en)*(sum(set.en$pi*!set.en$sub.y)/n1))
    Ft.w.en = 2*TPRt.w.en* PPVt.w.en/(TPRt.w.en + PPVt.w.en)
    en.w.IPW[[i]] = data.frame(TPRt.w.en, FPRt.w.en, PPVt.w.en, Ft.w.en, thres, nume.en.w, deno.en.w)
    
    ### Entropy - Ker
    e.proba.en.res2 = a+b*(set.en$sub.x %*%  tota.b)
    e.proba.en.res3 = a+b*(set.en$full.x %*%  tota.b)
    t.score = c(1 - 1 / (1 + exp(set.en$sub.x  %*%  tota.b)), 1 - 1 / (1 + exp(set.en$full.x %*%  tota.b)))
    r = c(rep(1,length(e.proba.en.res2)),rep(0,length(e.proba.en.res3)))
    bw <- npcdensbw(formula=factor(r)~c(t.score), bws = c(0, (t.n)^(-band)), bandwidth.compute= FALSE)
    data.eval.en <- data.frame(r=factor(rep(1, length(e.proba.en.res2))),t.score=1 - 1 / (1 + exp(set.en$sub.x  %*%  tota.b)))
    e.proba.en.res.ker = 1/fitted(npcdens(bws=bw, newdata = data.eval.en))
    
    deno.w.ker = sum(set.en$sub.y*e.proba.en.res.ker)/n1
    nume.w.ker = colSums(class.y.en*set.en$sub.y*e.proba.en.res.ker)/n1
    nume.w.f.ker =  colSums(e.proba.en.res.ker*class.y.en*!set.en$sub.y)/n1 
    TPR.w.en.ker = nume.w.ker/deno.w.ker
    FPR.w.en.ker = nume.w.f.ker/( 1 - deno.w.ker)
    PPV.w.en.ker = nume.w.ker/colMeans(class.y.true); PPV.w.en.ker = ifelse(PPV.w.en.ker>1 , 1, PPV.w.en.ker); PPV.w.en.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.en.ker = (1-FPR.w.en.ker) * (sum(e.proba.en.res.ker*!set.en$sub.y)/n1)/((1-TPR.w.en.ker)*deno.w.ker + (1-FPR.w.en.ker)*(sum(e.proba.en.res.ker*!set.en$sub.y)/n1))
    F.w.en.ker = 2*TPR.w.en.ker* PPV.w.en.ker/(TPR.w.en.ker + PPV.w.en.ker)
    en.w.ker[[i]] = data.frame(TPR.w.en.ker, FPR.w.en.ker, PPV.w.en.ker, F.w.en.ker, thres, nume.w.ker, deno.w.ker)
    
    est.score.Avar = c(set.en$sub.x %*% tota.b)
    pre.ind = which(FPR.w.en.ker <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarSEMIker.en.F[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, F.w.en.ker[pre.ind], deno.w.ker + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'SEMI')
    AvarSEMIker.en.ppv[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, PPV.w.en.ker[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'SEMI')
    
    
    ### Entropy - Spline
    reg = glm(r~1+ns(t.score,df=deg),family=binomial)
    e.proba.en.res.ker = 1/predict(reg,newdata=data.frame(t.score=1 - 1 / (1 + exp(set.en$sub.x  %*%  tota.b))),type="response")
    
    deno.w.ker = sum(set.en$sub.y*e.proba.en.res.ker)/n1
    nume.w.ker = colSums(class.y.en*set.en$sub.y*e.proba.en.res.ker)/n1
    nume.w.f.ker =  colSums(e.proba.en.res.ker*class.y.en*!set.en$sub.y)/n1 
    TPR.w.en.ker = nume.w.ker/deno.w.ker
    FPR.w.en.ker = nume.w.f.ker/( 1 - deno.w.ker)
    PPV.w.en.ker = nume.w.ker/colMeans(class.y.true); PPV.w.en.ker = ifelse(PPV.w.en.ker>1 , 1, PPV.w.en.ker); PPV.w.en.ker[colMeans(class.y.true) == 0] = 0
    NPV.w.en.ker = (1-FPR.w.en.ker) * (sum(e.proba.en.res.ker*!set.en$sub.y)/n1)/((1-TPR.w.en.ker)*deno.w.ker + (1-FPR.w.en.ker)*(sum(e.proba.en.res.ker*!set.en$sub.y)/n1))
    F.w.en.ker = 2*TPR.w.en.ker* PPV.w.en.ker/(TPR.w.en.ker + PPV.w.en.ker)
    en.w.spl[[i]] = data.frame(TPR.w.en.ker, FPR.w.en.ker, PPV.w.en.ker, F.w.en.ker, thres, nume.w.ker, deno.w.ker)
    
    est.score.Avar = c(set.en$sub.x %*% tota.b)
    pre.ind = which(FPR.w.en.ker <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarSEMIspl.en.F[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, F.w.en.ker[pre.ind], deno.w.ker + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'SEMI')
    AvarSEMIspl.en.ppv[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, PPV.w.en.ker[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'SEMI')
    
    est.score.Avar = c(set.en$sub.x %*% tota.b)
    pre.ind = which(FPRt.w.en <=A.thres)[1]; class.y.true.v = colMeans(e.proba.t>thres[pre.ind])
    AvarIPW.en.F[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, Ft.w.en[pre.ind], deno.en.w + class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'F', type1 = 'IPW')
    AvarIPW.en.ppv[i] = A.var(est.score.Avar, set.en$pi, set.en$sub.y, PPVt.w.en[pre.ind], class.y.true.v, class.y.true.v, n, thres[pre.ind], type = 'PPV', type1 = 'IPW')
    
    print(c(simple_prauc(Uni.IPW[[i]][,1], Uni.IPW[[i]][,2], Uni.IPW[[i]][,3]),
            simple_prauc(Uni.ker[[i]][,1], Uni.ker[[i]][,2] , Uni.ker[[i]][,3]),
            simple_prauc(Uni.spl[[i]][,1], Uni.spl[[i]][,2] , Uni.spl[[i]][,3])))
    print(c(simple_prauc(en.w.IPW[[i]][,1], en.w.IPW[[i]][,2], en.w.IPW[[i]][,3]),
            simple_prauc(en.w.ker[[i]][,1], en.w.ker[[i]][,2] , en.w.ker[[i]][,3]),
            simple_prauc(en.w.spl[[i]][,1], en.w.spl[[i]][,2] , en.w.spl[[i]][,3])))
    
  r.t = suppressMessages(roc(c(test.pool.y), c(e.proba.t)))
  t.measure[[i]] = as.matrix(coords(r.t, x = "all", input = "threshold", ret = c('tpr', "fpr",'ppv','npv','tp','fp','threshold')));
  auc.t[i] = suppressMessages(auc(test.pool.y, c(e.proba.t)))
  t.prauc1[i] = PRAUC(e.proba.t, test.pool.y)
  print(c( t.measure[[i]][which.max(t.measure[[i]][,7]>=A.thres),1] , 
           t.measure[[i]][which.max(t.measure[[i]][,7]>=A.thres),3] ,t.prauc1[i]))
  propt[i] = mean(test.pool.y)
  print(i)
} 
