source("PSweight_sub.R")

source("print_PSWeight_sub.R")
source("summary_PSWeight_sub.R")
source("print_PSweightsum_sub.R")
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

### example call - main effect glm with OW
# R option specified the number of bootstrap sample for SE calculation
p1<-PSweight_sub(ps.formula=ps.form_m,subgroup=subgroup,yname="Y",data=data,R=50,weight="overlap")

# can specify arbitary contrast statement
summary(p1,contrast = rbind(c(1,-1),c(0.5,-1)))

summary(p1,het = T)


# ATE - glm with IPW
p2 <- PSweight_sub(ps.formula=ps.form_m,subgroup=subgroup,yname="Y",data=data,R=50,weight="IPW")
p2
summary(p2,contrast = rbind(c(1,-1),c(0.5,-1)))

### ATT
p3 <- PSweight_sub(ps.formula=ps.form_m,subgroup=subgroup,yname="Y",trtgrp=0,data=data,R=50,weight="treated")
p3
summary(p3)

#entropy is also available
p4 <- PSweight_sub(ps.formula=ps.form_m,subgroup=subgroup,yname="Y",data=data,R=50,weight="entropy")
p4
summary(p4)

### can also include all variables from data in the ps model
p1<-PSweight_sub(ps.formula=Treatment ~.,subgroup=subgroup,yname="Y",data=data,R=50,weight="overlap")
summary(p1,contrast = rbind(c(1,-1),c(0.5,-1)))

### more flexible ps models
# pLASSO
p1<-PSweight_sub(ps.formula=ps,form_f,subgroup=subgroup,yname="Y",data=data,R=50, method='pLASSO',weight="overlap")
# can specify arbitary contrast statement
summary(p1,contrast = rbind(c(1,-1),c(0.5,-1)))
summary(p1,het = T)


# GBM: additional hyperparameter specification goes into the ps.control_sub
p1<-PSweight_sub(ps.formula=Treatment ~.,subgroup=subgroup,yname="Y",data=data,R=50, method='gbm',weight="overlap",ps.control_sub = list(n.trees=30))
summary(p1,contrast = rbind(c(1,-1),c(0.5,-1)))
summary(p1,het = T)

# RCS: the splined terms are included in the ps form
ps.form_rcs<-"Treatment~ X1+X2+X3+X4+X5+X6+X7+X8+rcs(X9)+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20"
p <- PSweight_sub(ps.formula = ps.form_rcs, subgroup = subgroup,yname = "Y", data=data,R=50, weight=c("overlap"), method = 'rcs')
summary(p,contrast = rbind(c(1,-1),c(0.5,-1)))
summary(p,het = T)


# RFs
p <- PSweight_sub(ps.formula = ps.form_m, subgroup = subgroup,yname = "Y", data=data, weight=c("overlap"), method = 'RFs')
summary(p,contrast = rbind(c(1,-1),c(0.5,-1)))
summary(p,het = T)


##################### user specified propensity score #############################################
### generate propensity
ps.estimate = p$propensity
p11<-PSweight_sub(ps.estimate =ps.estimate,zname='Treatment',subgroup=subgroup,yname="Y",data=data,R=50,weight="overlap")
summary(p11,contrast = rbind(c(1,-1),c(0.5,-1)))

### ATE
p2 <- PSweight_sub(ps.estimate =ps.estimate,zname='Treatment',subgroup=subgroup,yname="Y",data=data,R=50,weight="IPW")
p2
summary(p2,contrast = rbind(c(1,-1),c(0.5,-1)))

#a single row of propensity also supported
p2 <- PSweight_sub(ps.estimate =ps.estimate[,2],zname='Treatment',subgroup=subgroup,yname="Y",data=data,R=50,weight="IPW")
p2
summary(p2,contrast = rbind(c(1,-1),c(0.5,-1)))


### ATT
p3 <- PSweight_sub(ps.estimate =ps.estimate,zname='Treatment',subgroup=subgroup,yname="Y",trtgrp=0,data=data,R=50,weight="treated")
p3

summary(p3)

#entropy is also available
p4 <- PSweight_sub(ps.estimate =ps.estimate,zname='Treatment',subgroup=subgroup,yname="Y",data=data,R=50,weight="entropy")
p4
summary(p4)



