#load all the libraries required
library(Hotelling)
library(plyr)
library(mvtnorm)
library(nlme)
library(reshape2)
library(Matrix)
library(mvnormtest)

#load data
data <- read.csv("glass_data.csv",header=TRUE)
dim(data) #210 21 #1st 3 entries are the suspect and next 3 are crime scene rest background
#select only required eleents
glass.data <- data[,c(-4,-5,-6,-7,-9,-12,-19,-20)]
alpha = 0.1

head(glass.data)
head(data)

#We will keep only the following 10 elements for all the data calculations
# "K39" "Ti49" "Mn55" "Rb85" "Sr88" "Zr90" "Ba137" "La139" "Ce140" "Pb208"

#data partition
data.glass <- glass.data[,c(-1,-3)]
l.data <- log(data.glass[,-1],base=10)
colnames(l.data) <- paste("log", colnames(l.data), sep = "_") #modify column names
# log data is l.data
l.glass.data <- cbind(data.glass[,1],l.data)
colnames(l.glass.data) <- c('id', colnames(l.data))
head(l.glass.data)

#log data 
U <- l.glass.data[l.glass.data$id==1,]
S <- l.glass.data[l.glass.data$id==2,]
A <- l.glass.data[!(l.glass.data$id==1|l.glass.data$id==2),]

#function to check normality of each column of data
#normality of U
SW.norm <- function(data,alpha=0.05){
SW.statistic <- matrix(nrow=length(data), ncol=1)
SW.pvalue <- matrix(nrow=length(data), ncol=1)
SW.normality = vector("logical", length(data))
for (a in 1:length(data)){
  SW.statistic[a] <- shapiro.test(data[,a])$statistic
  SW.pvalue[a] <- shapiro.test(data[,a])$p.value  
  SW.normality[a] <- (SW.pvalue[a]>alpha)
}
result <- data.frame(colnames(data),SW.statistic,SW.pvalue,SW.normality)
colnames(result) <- c('Element','test statistic' , 'p-value','normality')
return(result)
}

#normality of U
SW.norm(glass.data[glass.data$fragment==1,-(1:3)],0.05) 
# everything normal except log_La139
SW.norm(U[,-1],0.05) 
# log scale everything except log_Sr88

#normality of S
SW.norm(glass.data[glass.data$fragment==2,-(1:3)],0.05) 
#everything normal
SW.norm(S[,-1],0.05)
#everything normal

#normality of A
SW.norm(glass.data[!(glass.data$fragment==1|glass.data$fragment==2),-(1:3)],0.05) 
#Only Pb208 normal
SW.norm(A[,-1],0.05) 
#log scale only 5 elements normal

#normality of whole data
SW.norm(glass.data[,-(1:3)],0.05) 
#Only Pb208 normal
SW.norm(l.glass.data[,-(1:3)],0.05) 
#log scale 5 elements  normal

######################################
#Approach 1: ASTM Two stage approach #
######################################

colnames(data.glass) <- c("id",colnames(data.glass[,-1]))

#we use the raw data for ASTM approach
U <- data.glass[data.glass$id==1,]
S <- data.glass[data.glass$id==2,]
A <- data.glass[!(data.glass$id==1|data.glass$id==2),]


Smean = colMeans(S[,-1])        # mean of each element
SDhat = apply(S[,-1], 2, sd)    # sd of each element
SD = cbind(SDhat, 0.03*Smean)   # (sd, minimum sd)                 
Imatch = cbind(Smean-4*apply(SD, 1, max), Smean+4*apply(SD, 1, max))   #match interval
Umean = colMeans(U[,-1])

P = length(Umean)
ASTMp = vector("logical", P)
for(i in 1:P){
  ASTMp[i] = (Imatch[i,1] <= Umean[i] & Umean[i] <= Imatch[i,2])
}
sum(ASTMp) == P
# TRUE means that the S is chemically indistinguishable from U
# Now we need to go to Stage 2 and find the RMP

### Let's start with the trace-anchored approach
N = length(unique(A$id))
ASTM_trace = matrix(nrow=N, ncol=P)
result_trace = vector("logical", N)
for(b in 1:N){
  Atemp = subset(A, A$id == unique(A$id)[b])[,-1]
  Amean = colMeans(Atemp)
  SDhatA = apply(Atemp, 2, sd)
  SDtemp = cbind(SDhatA, 0.03*Amean)
  Itemp = cbind(Amean-4*apply(SDtemp, 1, max), Amean+4*apply(SDtemp, 1, max))
  for(c in 1:P){
    ASTM_trace[b,c] = (Itemp[c,1] <= Umean[c] & Umean[c] <= Itemp[c,2])
  }
  result_trace[b] = (sum(ASTM_trace[b,]) == P)
  print(data.frame(item=unique(A$id)[b], match=result_trace[b]))
}
RMPtrace = mean(result_trace)
RMPtrace
# U would be chemically indistinguishable from 43.47% of the other sources


### Now let's try a source-anchored approach
ASTM_source = matrix(nrow=N, ncol=P)
result_source = vector("logical", N)
for(e in 1:N){
  Atemp = subset(A, A$id == unique(A$id)[e])[,-1]
  Amean = colMeans(Atemp)
  for(f in 1:P){
    ASTM_source[e,f] = (Imatch[f,1] <= Amean[e] & Amean[e] <= Imatch[f,2])
  }
  result_source[f] = (sum(ASTM_source[f,]) == P)
  print(data.frame(item=unique(A$id)[f], match=result_source[f]))
}
RMPsource = mean(result_source)
RMPsource
# 0% of the other sources would be chemically indistinguishable from S

## Last, let's try a general-match approach
ASTM_general = matrix(nrow=N*(N-1), ncol=P)
result_general = vector("logical", N*(N-1))
for(i in 1:N){
  Stemp = subset(A, A$id == unique(A$id)[i])[,-1]
  temp_Smean = colMeans(Stemp)
  temp_SDhat = apply(Stemp, 2, sd)
  SDtemp = cbind(temp_SDhat, 0.03*temp_Smean)
  ItempS = cbind(temp_Smean-4*apply(SDtemp, 1, max), temp_Smean+4*apply(SDtemp, 1, max))
  AnoS = subset(A, A$id != unique(A$id)[i])
  for(j in 1:(N-1)){
    Utemp = subset(AnoS, AnoS$id == unique(AnoS$id)[j])
    temp_Umean = colMeans(Utemp[,-1])
    count = j+(i-1)*(N-1)
    for(k in 1:P){
      ASTM_general[count,k] = (ItempS[k,1] <= temp_Umean[k] & temp_Umean[k] <= ItempS[k,2])
    }
    result_general[count] = (sum(ASTM_general[count,]) == P)
    print(data.frame(count, id1=unique(A$id)[i], id2=unique(AnoS$id)[j], match=result_general[count]))
  }
}
RMPgeneral = mean(result_general)
RMPgeneral
# 24.11% pairs of other sources would be chemically indistinguishable from eachother

#################################################
#Approach 2: Hotelling's T^2 Two stage approach #
#################################################
#We will consider dimension reduction for the hotelling test
#create principal components for dimension reduction
head(l.glass.data)
pc.data <- prcomp(l.glass.data[,-1], center = TRUE, scale. = TRUE)
print(pc.data)
summary(pc.data) #first 3 PC explain 96.38% of variation

pca_glass <- data.frame(
  id <- l.glass.data[,1],
  PC1 <- pc.data$x[,1],
  PC2 <- pc.data$x[,2],
  PC3 <- pc.data$x[,3]
)
colnames(pca_glass) <- c("id","PC1","PC2","PC3")
mshapiro.test(t(pca_glass)) 
# Shapiro-Wilk normality test
# data:  Z
# W = 0.94568, p-value = 9.857e-05
#reject multivariate normality

#test multivariate normality of 3 elements
# "id"        "log_K39"   "log_Ti49"  "log_Mn55"  "log_Rb85"  "log_Sr88" 
# "log_Zr90"  "log_Ba137" "log_La139" "log_Ce140" "log_Pb208"
mshapiro.test(t(l.glass.data[,c(8,9,11)]))
#we will consider "log_Ba137", "log_La139" and "log_Pb208" 
#for the Hotelling T2 approach
glass.data.mvn <- l.glass.data[,c(1,8,9,11)]
head(glass.data.mvn)

# seperate the suspect, crime scene and background data
U1 <- glass.data.mvn[glass.data.mvn$id==1,-1]
S1 <- glass.data.mvn[glass.data.mvn$id==2,-1]
A1 <- glass.data.mvn[!(glass.data.mvn$id==1|glass.data.mvn$id==2),]

#Stage 1
a <- hotelling.test(U1[,-1],S1[-1]) 
a$stats$statistic #test statistic = 15.47999
a$pval         #p-value 0.09304794
# at 5% level fail to reject equal means which 
# indicates there is an association in U and S

# Trace anchored probability
# Now we compare the remaining 23 windows WITH Eu
N <- length(unique(A1$id))
DU.hotelling.statistic <- matrix(nrow=N, ncol=1)
DU.hotelling.pvalue <- matrix(nrow=N, ncol=1)
DU.p = vector("logical", N)

for(j in 1:N){
  A1temp <- subset(A1, A1$id == unique(A1$id)[j])[,-1]
  temp <- hotelling.test(U1,A1temp)
  DU.hotelling.statistic[j] <- temp$stats$statistic
  DU.hotelling.pvalue[j] <- temp$pval
  DU.p[j] <- (DU.hotelling.pvalue[j]>alpha)
}
#trace anchored coincidence probability
trace.anchored.probability = sum(DU.p)/length(DU.p) #0.173913
#for alpha=0.1, we had association in stage 1 and small coincidence
#probability (17.39%) by trace anchored match approach which indicates 
#moderate evidence in favor of H_p


# Source anchored probability
# Now we compare the remaining 23 windows WITH Es
DS.hotelling.statistic <- matrix(nrow=N, ncol=1)
DS.hotelling.pvalue <- matrix(nrow=N, ncol=1)
DS.p = vector("logical", N)
for(j in 1:N){
  A2temp <- subset(A1, A1$id == unique(A1$id)[j])[,-1]
  temp <- hotelling.test(S1,A2temp)
  DS.hotelling.statistic[j] <- temp$stats$statistic
  DS.hotelling.pvalue[j] <- temp$pval
  DS.p[j] <- (DS.hotelling.pvalue[j]>alpha)
}
#source anchored coincidence probability
source.anchored.probability = sum(DS.p)/length(DS.p) #0.3043478

#for alpha=0.1, we had association in stage 1 and small coincidence
#probability (30.43%) by source anchored match approach which indicates 
#strong evidence in favor of H_p

#General match probability
# Now we compare the remaining 23 windows WITH Es
# for general match probability
# Total 23C2 = 253 combinations

GM.hotelling.statistic <- matrix(nrow=N, ncol=N)
GM.hotelling.pvalue <- matrix(nrow=N, ncol=N)
GM.p = matrix(FALSE, nrow=N, ncol=N)
for(i in 1:N){
  tempa <- subset(A1, A1$id == unique(A1$id)[i])[,-1]
  for(j in 1:N){
    tempb <- subset(A1, A1$id == unique(A$id)[j])[,-1]
    temp <- hotelling.test(tempa,tempb)
    GM.hotelling.statistic[i,j] <- temp$stats$statistic
    GM.hotelling.pvalue[i,j] <- temp$pval
    GM.p[i,j] <- (GM.hotelling.pvalue[i,j]>alpha)
  }
}
#general match coincidence probability
general.match.probability = (sum(GM.p)-N)/(N*(N-1)/2)  #0.229249
#for alpha=0.1, we had association in stage 1 and a moderately small coincidence
#probability (22.92%) by general match approach which indicates 
#strong evidence in favor of H_p

###################################################
#Approach 3: Likelihood Ratio on only 3 variables #
###################################################
#data partition
head(glass.data)
data.glass.lr <- glass.data[,c(2,3,10,11,13)]
#we will consider "log_Ba137", "log_La139" and "log_Pb208" 
l.data.lr <- log(data.glass.lr[,-c(1,2)],base=10)
#l.data.lr <- data.glass.lr[,-c(1,2)]
colnames(l.data.lr) <- paste("log", colnames(l.data.lr), sep = "_") #modify column names
l.glass.data.lr <- cbind(data.glass.lr[,c(1,2)],l.data.lr)
colnames(l.glass.data.lr) <- c('id','rep', colnames(l.data.lr))
head(l.glass.data.lr)
#Specific source NP LR
#partition the data
U <- l.glass.data.lr[l.glass.data.lr$id=='1',-2] #3 measurments
S <- l.glass.data.lr[l.glass.data.lr$id=='2',-2] #3 measurments
A <- l.glass.data.lr[!(l.glass.data.lr$id=='1'|l.glass.data.lr$id=='2'),-2]
#for u and s 
nu <- length(U$id) ; ns <- length(S$id) ; 
na <- length(unique(A$id)); nw <- length(A$id)/length(unique(A$id))
sum.s <- apply(S[,-1],2,sum) ; sum.u <- apply(U[,-1],2,sum)
theta.s <- (sum.u + sum.s)/(nu+ns);
sigma.s <- (((sum.u-theta.s)%*%t(sum.u-theta.s))+((sum.s-theta.s)%*%t(sum.s-theta.s)))/(nu+ns-1) 

#mean
glass_means <- ddply(l.glass.data.lr, c("id"), function(x) return(colMeans(x[,-(1:2)])))
glass_obs <- ddply(l.glass.data.lr, c("id","rep"), function(x) return(colMeans(x[,-(1:2)])))

set.seed(2323)
Eu <- data.frame(rep('u', nu), U[,-1])
colnames(Eu) = c('id', colnames(glass_means)[-1])

Es <- data.frame(rep('s', nu), S[,-1])
colnames(Es) = c('id', colnames(glass_means)[-1])

Ea <- data.frame(A)
colnames(Ea) = c('id', colnames(glass_means)[-1])

source('~/Reference code/BF.est.covs.R', chdir = TRUE)
source('~/Reference code/BF.make.sig.R', chdir = TRUE)

### Specific Source NPLR ###
mu_s_star = colMeans(rbind(Eu, Es)[,-1])
sig_s_star = cov(rbind(Eu, Es)[,-1])
temp_su = apply(rbind(Eu, Es)[,-1], 1, dmvnorm, mean=mu_s_star, sigma=sig_s_star, log=T)

mu_s = colMeans(Es[,-1])
sig_s = cov(Es[,-1])
temp_s = apply(Es[,-1], 1, dmvnorm, mean=mu_s, sigma=sig_s, log=T)

a_bar = ddply(Ea, c("id"), function(x) return(colMeans(x[,-1])))
mu_a = colMeans(a_bar[,-1])
mu_c_p = matrix(mu_a, 1, nw*3)
mu_c_p = as.vector(mu_c_p)
sig_a = s_star(Ea[,-1], Ea[,1])
sig_w = u_hat(Ea[,-1], Ea[,1])
sig_c_p = make.sig.c(sig_a, sig_w, nw)
long_a = matrix(t(Ea[,-1]), nrow=nw*3, ncol=length(unique(Ea$id)))
temp_a = apply(long_a, 2, dmvnorm, mean=mu_c_p, sigma=sig_c_p, log=T)

# a_bar_star = rbind(a_bar[,-1], colMeans(Eu[,-1]))
# mu_a_star = colMeans(a_bar_star)
# mu_c_d = matrix(mu_a_star, 1, nw*10)
# mu_c_d = as.vector(mu_c_d)
a_bar_star = rbind(a_bar[,-1], colMeans(Eu[,-1]))
mu_a_star = colMeans(a_bar_star)
# mu_c_d = matrix(mu_a_star, 1, nw*10)
# mu_c_d = as.vector(mu_c_d)
mu_c_d1 = matrix(mu_a_star, 1, nw*3)
mu_c_d1 = as.vector(mu_c_d1)
mu_c_d2 = matrix(mu_a_star, 1, nu*3)
mu_c_d2 = as.vector(mu_c_d2)
sig_a_star = s_star(rbind(Ea,Eu)[,-1], rbind(Ea, Eu)[,1])
sig_w_star = u_hat(rbind(Ea,Eu)[,-1], rbind(Ea, Eu)[,1])
#sig_a_star1 = s_star(Ea[,-1], Ea[,1])
#sig_w_star1 = u_hat(Ea[,-1], Ea[,1])
#sig_a_star2 = s_star(Eu[,-1], Eu[,1])
#sig_w_star2 = u_hat(Eu[,-1], Eu[,1])

sig_c_d1 = make.sig.c(sig_a_star, sig_w_star, nw)
sig_c_d2 = make.sig.c(sig_a_star, sig_w_star, nu)
#sig_c_d1 = make.sig.c(sig_a_star1, sig_w_star1, nw)
#sig_c_d2 = make.sig.c(sig_a_star2, sig_w_star2, nu)
###########################################################
#long_au = matrix(t(rbind(Ea,Eu)[,-1]), nrow=nw*3, ncol=length(unique(Ea$id))+1)
long_au1 = matrix(t(Ea[,-1]), nrow=nw*3, ncol=length(unique(Ea$id)))
long_au2 = matrix(t(Eu[,-1]), nrow=nu*3, ncol=1)
#temp_au = apply(long_au, 2, dmvnorm, mean=mu_c_d, sigma=sig_c_d, log=T)
temp_au1 = apply(long_au1, 2, dmvnorm, mean=mu_c_d1, sigma=sig_c_d1, log=T)
temp_au2 = apply(long_au2, 2, dmvnorm, mean=mu_c_d2, sigma=sig_c_d2, log=T)

numerSS = sum(temp_su) + sum(temp_a)
denomSS = sum(temp_s) + sum(temp_au1)+sum(temp_au2)
logNPLRss = numerSS - denomSS
exp(logNPLRss) #0.06561043 where logNPLRss = -2.724021

# LR is less than 1 which indicates more likely to observe
# the evidence if Hd is true than if Hp is true


################################################
#Approach 4: Compute Bayes Factors - using ABC #
################################################
library(mvtnorm)
library(MCMCglmm)
library(Matrix)
library(pdist)
background <- read.csv("background.csv",header=TRUE)
names(background)
BD <- background[,c(-3,-4,-5,-6,-8,-11,-18,-19)]
BD2 <- log(BD[,-c(1:2)],base=10)
colnames(BD2) <- paste("log", colnames(BD2), sep = "_")

#background data
bg.data <- cbind(background[,2],BD2)
colnames(bg.data) <- c('id', colnames(BD2))
head(bg.data)

# The data
Eu <- l.glass.data[l.glass.data$id==1,]
Es <- l.glass.data[l.glass.data$id==2,]
Ea <- l.glass.data[(l.glass.data$id!=1|l.glass.data$id!=2),]
D <- bg.data

################################################################################
#source files
source("~/Reference code/BF.est.MC.R")
source("~/Reference code/BF.approx.MCSE.R")
source("~/Reference code/BF.make.sig.R")
source("~/Reference code/BF.est.covs.R")

source("~/Reference code/BF.StdMean.SSsim.f1.R")
source("~/Reference code/BF.HarmMean.SSsim.f1.R")
source("~/Reference code/BF.StdMean.CSsim.f1.R")
source("~/Reference code/BF.HarmMean.CSsim.f1.R")

source("~/Reference code/BF.StdMean.SSsim.f2.R")
source("~/Reference code/BF.HarmMean.SSsim.f2.R")
source("~/Reference code/BF.StdMean.CSsim.f2.R")
source("~/Reference code/BF.HarmMean.CSsim.f2.R")
################################################################################

prior.mean=colMeans(D[,-1])
prior.within.cov=u_hat(D[,-1], D[,1])
prior.between.cov=s_star(D[,-1], D[,1])
prior_list55=list(mean=prior.mean, sigb=prior.between.cov, dofb=55, 
                  sige=prior.within.cov, dofe=55, K=10)

D_means = ddply(D, c('id'), function(x){ return(colMeans(x[,-1])) })
summary(dist(D_means[,-1]))
tol = 0.4

# for the specific source
BF.ABC.SS <-function(dat_u, dat_s, dat_a, prior_theta, prior_Mp, epsilon, niters=30000, burn_in=1000, nthin=15){
  
  # MAKE SURE THAT dat_u, dat_s, and dat_s ALL HAVE THE SAME NUMBER OF COLUMNS!!
  ncol_u = dim(dat_u)[2]
  ncol_s = dim(dat_s)[2]
  ncol_a = dim(dat_a)[2]
  if(ncol_u != ncol_s){print("Data Dimension Error")}
  if(ncol_u != ncol_a){print("Data Dimension Error")}
  num_vars=dim(dat_u)[2]-1
  
  ### Sample m from prior_model
  Np = rbinom(1, niters, prior_Mp)
  Nd = niters - Np
  
  ###
  # PRIOR - prosecution model
  # B: prior for mu_s ~ MVN(mu, V_b)
  # R: prior for sig_s ~ InvWish(V_e, nu_e)
  pr<-list(B=list(mu=prior_theta$mean, V=prior_theta$sigb),
           R=list(V=prior_theta$sige, nu=prior_theta$dofe))
  
  # Sample values for the parameters mu_s and sig_s from the posterior distribution - pi(theta|e_s)
  # Trait and units are internal functions for doing the multivariate models.
  theta_s<-MCMCglmm(as.formula(paste("cbind(",paste(colnames(dat_s)[-1],collapse=","),") ~ trait-1", sep="")),
                    rcov = ~us(trait):units, data = dat_s[,-1], family = rep("gaussian", num_vars),
                    verbose = FALSE, prior=pr, nitt=nthin*Np + burn_in, burnin=burn_in, thin=nthin)
  
  # Grab the samples for mu_s and sig_s from the Gibbs sampler output
  mu_s=theta_s$Sol
  
  sig_s<-array(NA, c(num_vars, num_vars, Np))
  for (i in 1:Np){
    sig_s[,,i]<-matrix(theta_s$VCV[i,], num_vars, num_vars)
  }
  
  ###
  # PRIOR - defense model
  # B: prior for mu_a ~ MVN(mu, K*V_b)
  # R: prior for sig_w ~ InvWish(V_e, nu_e)
  # G: prior for sig_a ~ InvWish(V_b, nu_b)		
  pr<-list(B=list(mu=prior_theta$mean, V=prior_theta$K*prior_theta$sigb),
           R=list(V=prior_theta$sige, nu=prior_theta$dofe),
           G=list(G1=list(V=prior_theta$sigb, nu=prior_theta$dofb)))
  
  # Sample values for the parameters mu_a_prior, sig_a_prior, sig_w_prior 
  # from the posterior distribution - pi(theta|e_a)
  # Trait and units are internal functions for doing the multivariate models.			
  theta_a<- MCMCglmm(as.formula(paste("cbind(",paste(colnames(dat_a)[-1],collapse=","),") ~ trait-1", sep="")),
                     random =as.formula(paste("~us(trait):",colnames(dat_a)[1], sep="")),
                     rcov = ~us(trait):units, data = dat_a, family = rep("gaussian", num_vars),
                     verbose = FALSE, prior=pr, nitt=nthin*Nd + burn_in, burnin=burn_in, thin=nthin)
  
  # Grab the samples for mu_a, sig_a, sig_w from the Gibbs sampler output
  mu_a=theta_a$Sol
  
  sig_a<-array(NA, c(num_vars, num_vars, Nd))
  for (i in 1:Nd){
    sig_a[,,i]<-matrix(theta_a$VCV[i,1:(num_vars* num_vars)], num_vars, num_vars)
  }
  
  sig_w<-array(NA, c(num_vars, num_vars, Nd))
  for (i in 1:Nd){
    sig_w[,,i]<-matrix(theta_a$VCV[i,(num_vars* num_vars+1):(2*num_vars* num_vars)], num_vars, num_vars)
  }
  
  ### ABC with Euclidean distance between mean vectors and equal model priors
  accept_Mp = vector("logical", Np)
  for(i in 1:Np){
    u_star = rmvnorm(nrow(dat_u), mean=mu_s[i,], sigma=sig_s[,,i])
    accept_Mp[i] = dist(rbind(colMeans(u_star), colMeans(dat_u[,-1]))) < epsilon
  }
  accept_Md = vector("logical", Nd)
  for(i in 1:Nd){
    A_star = rmvnorm(1, mean=mu_a[i,], sigma=sig_a[,,i])
    u_star = rmvnorm(nrow(dat_u), mean=A_star, sigma=sig_w[,,i])
    accept_Md[i] = dist(rbind(colMeans(u_star), colMeans(dat_u[,-1]))) < epsilon
  }
  
  S = sum(accept_Mp) + sum(accept_Md)
  BF_num = sum(accept_Mp)/prior_Mp
  BF_den = sum(accept_Md)/(1-prior_Mp)
  BF_abc = c(S, BF_num, BF_den, BF_num/BF_den)
  
  ###
  parameters = list(mu_s, sig_s, mu_a, sig_a, sig_w)
  names(parameters) = c("mu.s", "Sig.s", "mu.a", "Sig.a", "Sig.w")
  return_list = list(BF_abc, parameters)
  names(return_list) = c("bf", "theta")
  return(return_list)
}

#ABC Bayes factor 
set.seed(23)
BF.ABC.SS(Eu, Es, Ea, prior_list55, 0.2, tol)$bf
#14019.000000 31030.000000 12128.888889     2.558355
BF.ABC.SS(Eu, Es, Ea, prior_list55, 0.5, tol)$bf
#20998.000000 29760.000000 12236.000000     2.432167
BF.ABC.SS(Eu, Es, Ea, prior_list55, 0.8, tol)$bf
#28218.00000 29994.44444 12230.00000     2.45253