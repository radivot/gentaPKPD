#This script does some statistics using the results of multiRegimenPKPD.R (see also comments in PKPD_6.R)
#Outputs are the rho and p values of the paper and Fig. 5.
rm(list=ls())
inn<-read.csv("statisticsinput.csv",header=F)#statistics.csv gets its numbers from sumOutput.csv
rrho=(with(inn,cor.test(V1,V2,method = "spearman"))) # TR: this probably suffices for this section
print(rrho)
# Spearman's rank correlation rho
# data:  V1 and V2
# S = 378, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates: # rho=0.9307185 

two<-as.vector(100*(inn[,2]))#each entry is sum of 2 probabilities(%) and therefore in [0,200] 
predictor<-1:32 #lexicograhic order, merely ordinal, not quantitative.
plot(predictor,two, xlim=c(0,32), pch=19,bty="l",ann=0)# Fig. 5 of our paper