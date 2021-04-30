

# source functions
# set working directory
source("SumStat_sub.R")
source("print_Sumstat_sub.R")
source("Plot_SumStat_sub.R")
source('PSmethod_sub.R')
#########################################
######## Illustration 1 #################

data <- read.csv("sample_data.csv")[,-1]
# pre-specify the confounder names: X1-X20
cov <- paste0('X',1:20)


# pre-specify subgroups by column index or column names
subg <- data[,c(2,3,4,20,21)]
subgroup<-paste0('X',c(1,2,3, 19,20))


# specify a main effect propensity score model
ps.form_m <- paste("Treatment~",paste0("X",1:20,collapse = "+"))
# or any arbitrary ps model
ps.form1 <- paste(ps.form_m,"+X1*X7+X1*X8+X2*X7")
# or a full interaction ps model
ps.form_f = ps.form_m
for(i in 1:length(subgroup)){
  for (j in 1:length(cov)){
    ps.form_f=paste(ps.form_f, paste(subgroup[i],"*",cov[j]), sep=" + ")
  }
}

# can also provide the covariate matrix  
covM<-as.data.frame(model.matrix(formula(ps.form_m),data))[,-1]

##################### user specified propensity score model #############################################

### example call - main effect glm with OW
p <- SumStat_sub(ps.formula = ps.form_m, subgroup = subgroup, data=data, weight=c("overlap"))
p

# default x axis print all main effect variables, excluding subgroup variables
plot(p)
# the base= T option prints unadjusted ASMDs
plot(p, base = T)
# the plotsub = T option prints subgroup variables on x axis as well
plot(p, plotsub = T)

# can specify certain variables by index or variable names
plot(p,varlist = c(1,2,3))
plot(p,varlist = c("X1","X2","X7"))
plot(p,varlist = c("X1","X2","X7"),base = T)


# ATE - glm with IPW

p <- SumStat_sub(ps.formula = ps.form2, subgroup = subgroup, data=data, weight=c("IPW"))
p
plot(p)

# ATT - glm with ATT weights
p <- SumStat_sub(ps.formula = ps.form2, subgroup = subgroup, data=data, weight=c("treated"))
p
plot(p)

# ATEN - glm with entropy weights
p <- SumStat_sub(ps.formula = ps.form2, subgroup = subgroup, data=data, weight="entropy", trtgrp = 0)
p
plot(p)


### can also include all variables from data in the ps model
p <- SumStat_sub(ps.formula = Treatment~., subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"))
p

### more flexible ps models
# pLASSO
p <- SumStat_sub(ps.formula = ps.form_f, subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"),method = 'pLASSO')
p
plot(p,plotsub = T)

# GBM: additional hyperparameter specification goes into the ps.control_sub
p <- SumStat_sub(ps.formula = Treatment~., subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"),method = 'gbm',ps.control_sub = list(n.trees=30))
p

# RCS: the splined terms are included in the ps form
ps.form_rcs<-"Treatment~ X1+X2+X3+X4+X5+X6+X7+X8+rcs(X9)+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20"
p <- SumStat_sub(ps.formula = ps.form_rcs, subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"), method = 'rcs')
p
# this plots the ASMD of the specified formula
plot(p)

# specify covariate matrix covM can display only covariates
p <- SumStat_sub(ps.formula = ps.form_rcs, subgroup = subgroup,yname = "Y", data=data, covM =covM,  weight=c("overlap"), method = 'rcs')
plot(p)


# RFs
p <- SumStat_sub(ps.formula = ps.form_m, subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"),covM=covM, method = 'RFs')
p
plot(p)
##################### user specified propensity score #############################################
### generate propensity
p0 <- SumStat_sub(ps.formula = ps.form_m, subgroup = subgroup, data=data, weight=c("overlap"))
propfull<-p0$propensity
proprow<-propfull[,1]
trtgrp<-1 
# default print all variables in data
plot(p)

### can specify covariate matrix to display by using the covM option
p0 <- SumStat_sub(ps.formula = ps.form_m, subgroup = subgroup, data=data,covM = covM, weight=c("overlap"))
plot(p0)

# ATO
p <- SumStat_sub(zname = "Treatment", ps.estimate = propfull,trtgrp=trtgrp, subgroup = subgroup, data=data, covM=covM, weight=c("overlap"))
p
plot(p)
plot(p,base = T)


# ATE
p <- SumStat_sub(zname = "Treatment",ps.estimate = propfull,trtgrp=trtgrp, subgroup = subgroup,  covM=covM,data=data, weight=c("IPW"))
plot(p)


# ATT
p <- SumStat_sub(zname = "Treatment", ps.estimate = proprow,trtgrp=0,subgroup = subgroup, data=data, covM=covM, weight=c("treated"))
plot(p)


