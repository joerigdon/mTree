# from https://github.com/joerigdon/HTE/blob/master/GenData1A.R
library(Hmisc)

rmulti = function(prob, n) {
  as.numeric(rMultinom(probs=matrix(prob, 1, length(prob)), m=n))
}


##Set up data where B~bern(0.3), eGFR~norm(60, 20), age~norm(60,10), sex~bern(0.5)
true = data.frame(age=rnorm(200, 65, 5), stage4=factor(rbinom(200, 1, 0.98)), site=factor(rmulti(c(0.5, 0.17, 0.33), 200)), prevTrt=factor(rmulti(c(0.5, 0.2, 0.2, 0.1), 200)), ecog=factor(rbinom(200, 1, 0.55)), diseaseFree=factor(rbinom(200, 1, 0.35)))

true$delta = rmulti(c(0.5, 0.2, 0.3), 200)-2
mean(true$delta)

true$Y1 = NA
true$Y0 = NA

true$Y1[true$delta==-1] = 0
true$Y0[true$delta==-1] = 1

true$Y1[true$delta==1] = 1
true$Y0[true$delta==1] = 0

true$Y1[true$delta==0] = 1
true$Y0[true$delta==0] = 1

sum((true$Y1-true$Y0)!=true$delta) #they all equal delta
mean(true$Y1-true$Y0) #-0.15433 as specified above for ATE
mean(true$delta)
mean(true$Y1)
mean(true$Y0)


##Add in rowname
true$UNIQID = paste("A", rownames(true), sep="")

##Change types where necessary
true$delta = factor(true$delta)
true$Y1 = factor(true$Y1)
true$Y0 = factor(true$Y0)

# Take random sample of potential outcomes to get experiment,
set.seed(12)
ind.trt = sample(200, 100)
true$Z = 0
true$Z[ind.trt] = 1
true$Y = NA
true$Y[true$Z==1] = true$Y1[true$Z==1]
true$Y[true$Z==0] = true$Y0[true$Z==0]
true$Y = true$Y-1

example_trial_data <- true
##Save
usethis::use_data(example_trial_data, overwrite = TRUE)




