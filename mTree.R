##Load packages/dependencies
library(optmatch)
library(party)
library(MASS)
library(Gmisc, verbose=FALSE)
library(Hmisc)
library(RI2by2)
library(LTRCtrees)
library(survival)


##0. Helper functions
##Function to record balance
asd = function(v, Z, ind.cat) {
 v = as.numeric(as.character(v))
 trt = v[Z==1]
 ctrl = v[Z==0]
 m.trt = mean(trt, na.rm = TRUE)
 m.ctrl = mean(ctrl, na.rm = TRUE)
 sd2 = sqrt(var(trt, na.rm = TRUE)/2 + var(ctrl, na.rm = TRUE)/2)
 asd = abs(m.trt - m.ctrl)/sd2
 if (ind.cat==1) {
  asd = double()
  tt = as.numeric(names(table(v)))
  for (i in 1:length(tt)) {
   m.trt = sum(trt==tt[i])/length(trt)
   m.ctrl = sum(ctrl==tt[i])/length(ctrl)
   sd2 = sqrt(m.trt * (1 - m.trt)/2 + m.ctrl * (1 - m.ctrl)/2)
   asd = c(asd, abs(m.trt - m.ctrl)/sd2)
  }
 }
asd
}


##Function to make tables of balances
mktabASD = function(data, var.names, ind.cat, group.name, miss, digit) {
 n = names(data)
 cols = data[, which(n==group.name)]
 r = table(data[, which(n==group.name)])
 j = c(apply(r, 1, function(x) paste("n=",x,sep="")))
 for (i in 1:length(var.names)) {
  if (ind.cat[i]==0) {
   outc = data[, which(n==var.names[i])]
   label(outc) = var.names[i]
   dd2 = getDescriptionStatsBy(outc, cols, html=TRUE, useNA=miss, statistics=FALSE, add_total_col=FALSE, continuous_fn=describeMean, digits=digit)
   dd2 = cbind(dd2, round(asd(outc, cols, 0), 2))
   rownames(dd2)[1] = var.names[i]
  } else if (ind.cat[i]==1) {
     outc = data[, which(n==var.names[i])]
     label(outc) = var.names[i]
     dd1 = getDescriptionStatsBy(factor(outc), cols, html=TRUE, useNA=miss, statistics=FALSE, add_total_col=FALSE, digits=digit)
     dd2 = rbind(rep("", dim(dd1)[2]), dd1)
     dd2 = cbind(dd2, c("", round(asd(outc, cols, 1), 2)))
     rownames(dd2)[1] = var.names[i]
    }
    j = rbind(j,dd2)
  }
#Clean up a few things (orderings from categorical variables, etc.)
 k = t(apply(j, 1, function(x) gsub("&plusmn;", "\\±", x)))
 k2 = t(apply(k, 1, function(x) gsub("&lt; ", "<", x)))
 rownames(k2)[rownames(k2)=="j"] = ""
 rownames(k2) = lapply(rownames(k2),function(x) gsub("a[.]", " ", x))
 rownames(k2) = lapply(rownames(k2),function(x) gsub("b[.]", " ", x))
 rownames(k2) = lapply(rownames(k2),function(x) gsub("c[.]", " ", x))
 rownames(k2) = lapply(rownames(k2),function(x) gsub("d[.]", " ", x))
 rownames(k2) = lapply(rownames(k2),function(x) gsub("e[.]", " ", x))
 rownames(k2) = lapply(rownames(k2),function(x) gsub("f[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("a[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("b[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("c[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("d[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("e[.]", " ", x))
 colnames(k2) = lapply(colnames(k2),function(x) gsub("f[.]", " ", x))
#Remove uninformative missing values if missing="ifany"
 if (miss=="always") {
  n0 = apply(k2, 1, function(x) sum(x=="0 (0%)" | x=="0 (0.0%)"))
  rmv = which(rownames(k2)=="Missing" & n0==2)
  if (length(rmv)>0) {
   k2 = k2[-as.numeric(rmv), ]
  }
 }
 if (dim(k2)[2]>2) {
  k2[1, 3] = ""
  colnames(k2)[3] = "ASD"
 }
k2
}


##Function to make rank based Mahalanobis distance matrix
smahal = function(X){
 X = as.matrix(X)
 n = dim(X)[1]
 k = dim(X)[2]
 for (j in 1:k) {
  X[, j]=rank(X[, j]) #compute on ranks
 }
 cv = cov(X)
 vuntied = var(1:n)
 rat = sqrt(vuntied/diag(cv))
 cv = diag(rat) %*% cv %*% diag(rat)
 out = matrix(NA, n, n)
 icov = ginv(cv)
 for (i in 1:n) out[i, ] = mahalanobis(X, X[i, ],icov, inverted=T)
 out
}


##Function to create formulas as needed (do we need?)
formulize = function(outc, covs) {
 string = paste( c(outc, paste(covs, collapse=" + ") ), collapse=" ~ " )
 return(as.formula(string))
}


##Function to obtain estimates for binary setting
getEst_b = function(e1) {
 e1$Y = as.numeric(as.character(e1$Y))   
 p1 = mean(e1$Y[e1$Z==1])
 p0 = mean(e1$Y[e1$Z==0])
 ee1 = matrix(table(-e1$Z, -e1$Y), 2, 2, byrow=FALSE)
 r2 = Perm.CI.RLH(ee1, 0.05, total_tests=1000)
 data.frame(x1=ee1[1, 1], n1=sum(ee1[1, ]), p1=p1, x0=ee1[2, 1], n0=sum(ee1[2, ]), p0=p0, diff=r2$tau.hat, lower=r2$RLH[1], upper=r2$RLH[2])
}


##Function to obtain estimates for continuous setting
getEst_c = function(e1) {
 j = wilcox.test(e1$Y[e1$Z==1], e1$Y[e1$Z==0], conf.int=TRUE)
 data.frame(diff=j$estimate, lower=j$conf.int[1], upper=j$conf.int[2])
}


##Function to obtain (additive) estimates for survival setting
##Add multiplicative?
getEst_s = function(e1, l, u, l2=-10, u2=3) {
 est = seq(l, u, (u-l)/1000)
 pval = double()
 for (i in 1:length(est)) {
  dta2 = e1
  dta2$Y[dta2$Z==0] = dta2$Y[dta2$Z==0]+est[i]
  j2 = survdiff(Surv(Y, event) ~ Z, data=dta2, rho=0) #logrank
  pval = c(pval, j2$pvalue)
  }
 add = data.frame(diff=est[which.max(pval)], lower=min(est[pval>0.05]), upper=max(est[pval>0.05]))
 est2 = exp(seq(l2, u2, (u2-l2)/1000))
 pval2 = double()
 for (i in 1:length(est2)) {
  dta2 = e1
  dta2$Y[dta2$Z==0] = dta2$Y[dta2$Z==0]*est2[i]
  j2 = survdiff(Surv(Y, event) ~ Z, data=dta2, rho=0) #logrank
  pval2 = c(pval2, j2$pvalue)
  }
 mult = data.frame(diff=est2[which.max(pval2)], lower=min(est2[pval2>0.05]), upper=max(est2[pval2>0.05]))
 output.all = list(add=add, mult=mult)
 return(output.all)
}


##1. Matching step
##Function to create pair-match
mt_match = function(dta, id, trt, outc, covs, ind.cat, censVar=NA) {
 dta2 = dta[order(dta[, names(dta)==trt], decreasing=TRUE), ]
 trt2 = dta2[, names(dta2)==trt] #get trt variable
 id2 = dta2[, names(dta2)==id] #get ID variable
 dta2$tID = paste(trt2, id2, sep="_") #ID for match
#Convert factors to numeric where necessary
 dta2[, which(ind.cat==1)] = sapply(dta2[, which(ind.cat==1)], function(x) as.numeric(as.character(x)))
 dta3 = dta2[, names(dta2) %in% covs] #data for match
 dist = smahal(X=dta3) #Need numeric for match; could parallelize here
 rownames(dist) = dta2$tID
 colnames(dist) = dta2$tID
 dist2 = dist[substr(rownames(dist), 1, 1)==1, substr(colnames(dist), 1, 1)==0] #get distance matrix of treated (rows) by control (columns); important to get weights from here
 if (dim(dist2)[1]>dim(dist2)[2]) {
  dist2 = t(dist2)
 } #transpose if more treated than controls
 options("optmatch_max_problem_size" = Inf) #accommodate larger datasets
 pm = pairmatch(dist2) #execute pair-match; parallelize here
 org = data.frame(ID=names(pm), match=pm)
#Make into treated-control ID data frame
 org2 = org[order(org$match, org$ID), ] #stack on top of each other
 org3 = org2[!is.na(org2$match), ]
 nmatch = dim(org3)[1]/2
 jj = data.frame(con=rep(NA, nmatch), trt=rep(NA, nmatch), match=rep(NA, nmatch))
 for (i in 1:nmatch) {
  jj$con[i] = as.character(org3$ID[2*i-1])
  jj$trt[i] = as.character(org3$ID[2*i])
  jj$match[i] = as.character(org3$match[2*i])
 }
#Add weight to jj
 jj$wt = NA
 for (i in 1:dim(jj)[1]) {
  jj$wt[i] = dist2[rownames(dist2)==jj$trt[i], colnames(dist2)==jj$con[i]]
 }  
#Make delta data set
 t0 = dta2[match(jj$trt, dta2$tID), ] #match preserves ordering
 t0$match = jj$match[jj$trt==t0$tID]
 rownames(t0) = t0$match
 t1 = t0[, names(t0) %in% covs]
 c0 = dta2[match(jj$con, dta2$tID), ]
 c0$match = jj$match[jj$con==c0$tID]
 rownames(c0) = c0$match
 c1 = c0[, names(c0) %in% covs]
 aa = list(t1, c1)
 dd = data.frame(delta=t0[, names(t0)==outc]-c0[, names(c0)==outc], Reduce(`+`, aa) / length(aa)) #create delta; average covariates
#Fix those with censored outcomes
 if (!is.na(censVar)) {
  dd = dd[, names(dd)!="delta"]
  dd$match = rownames(dd)
  dd$timeT = NA
  dd$eT = NA
  dd$timeC = NA
  dd$eC = NA
  for (i in 1:dim(dd)[1]) {
   dd$timeT[i] = t0[which(t0$match==dd$match[i]), names(t0)==outc]
   dd$eT[i] = t0[which(t0$match==dd$match[i]), names(t0)==censVar]   
   dd$timeC[i] = c0[which(c0$match==dd$match[i]), names(c0)==outc]
   dd$eC[i] = c0[which(c0$match==dd$match[i]), names(c0)==censVar]   
  }
  dd$time1 = -Inf
  dd$time2 = Inf
  for (i in 1:dim(dd)[1]) {
   if (dd$eT[i]==1 & dd$eC[i]==1) {
    dd$time1[i]=dd$timeT[i]-dd$timeC[i]
    dd$time2[i]=dd$timeT[i]-dd$timeC[i]
   } else if (dd$eT[i]==1 & dd$eC[i]==0) {
      dd$time2[i]=dd$timeT[i]-dd$timeC[i]
     } else if (dd$eT[i]==0 & dd$eC[i]==1) {
        dd$time1[i]=dd$timeT[i]-dd$timeC[i]
       }
  }
 dd = dd[, match(c("time1", "time2", covs), names(dd))] 
 }
#Merge the weights onto dd
 dd$match = rownames(dd)
 dd2 = merge(dd, jj[, names(jj) %in% c("match", "wt")], by="match", all.x=TRUE)
 rownames(dd2) = dd2$match 
 dd3 = dd2[, names(dd2)!="match"]
#Remove the categorical variables where trt != con
 covs2 = covs[which(ind.cat==1)]
 t1cat = t1[, names(t1) %in% covs2]
 c1cat = c1[, names(c1) %in% covs2]
 elim = double() #want to equal 0 to keep
 if (length(covs2)==1) {
  elim = as.numeric(t1cat!=c1cat)
 } else if (length(covs2)>1) {
    elim = apply(t1cat!=c1cat, 1, sum)
   } 
#Create data frame, making sure categorical variables are factors
 dd4 = dd3[elim==0, ]
 for (i in 1:length(covs2)) {
   dd4[, names(dd4)==covs2[i]] = as.factor(dd4[, names(dd4)==covs2[i]])
 } 
 dd4
}


##2. Tree step
##Function to obtain decision tree
mt_tree = function(dta, covs, surv=FALSE, l=-Inf, r=Inf, wtVar=NA, G=4) {
 dta$weight = ifelse(is.na(wtVar), NA, (G+1)-as.numeric(cut2(dta[, names(dta)==wtVar], g=G)))
 if (surv==FALSE) {       
  if (dim(table(dta$delta))<=3) {dta$delta = as.factor(dta$delta)}
  if (is.na(wtVar)) {
   fit = ctree(formulize(outc="delta", covs=covs), data=dta)
  } else if (!is.na(wtVar)) {
     dta$weight = (G+1)-as.numeric(cut2(dta[, names(dta)==wtVar], g=G))
     fit = ctree(formulize(outc="delta", covs=covs), data=dta, weights=dta$weight)
    }  
 } else if (surv==TRUE) {
    dta$time1[dta$time1==-Inf] = l
    dta$time2[dta$time2==Inf] = r
    fit = LTRCIT(formulize(outc="Surv(time1, time2, type='interval2')", covs), data=dta)
   }  
fit 
}


##3. Estimation/statistics step
##Function to get balance stats and effect estimates within subgroups
mt_stats = function(dta, trt, outc, tree, covs, ind.cat, censVar=NA) {
#Make sure variable types are appropriate
 FAC = covs[which(ind.cat==1)]
 for (i in 1:length(FAC)) {
  dta[, names(dta)==FAC[i]] = as.factor(dta[, names(dta)==FAC[i]])
 }
 NUM = covs[which(ind.cat==0)]
 for (i in 1:length(NUM)) {
  dta[, names(dta)==NUM[i]] = as.numeric(dta[, names(dta)==NUM[i]])
 }    
#Get subgroup label for original observations
 dta$predNode = predict(tree, newdata=dta, type="node")    
#Look at balance within discovered subgroups
 dta$Z = dta[, names(dta)==trt]
 tabs = by(dta, dta$predNode, function(x) mktabASD(x, var.names=covs, ind.cat=ind.cat, group.name="Z", miss='always', digit=1))
#Look at estimates within discovered subgroups
 dta$Y = dta[, names(dta)==outc]
 if (length(unique(dta$Y))<=2 & is.na(censVar)) {   
  est = by(dta, dta$predNode, function(x) getEst_b(x))
 } else if (length(unique(dta$Y))>2 & is.na(censVar)) {   
  est = by(dta, dta$predNode, function(x) getEst_c(x))
 } else {
    dta$event = dta[, names(dta)==censVar]
    L = -max(dta$Y)
    U = max(dta$Y)
    est = by(dta, dta$predNode, function(x) getEst_s(x, L, U))
   }  
 output.all = list(tabs=tabs, est=est)
 return(output.all)  
}


