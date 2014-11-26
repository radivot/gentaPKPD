#This script depends on PKPD_6.R, whose comments apply also to this script and should be skimmed first. 
#Input is a list of many dosing regimens, not just one regimen.
#Output is a data frame, for all regimens input and for both GA=25 and GA=40. For figure output use PKPD_6.R
rm(list=ls())
source("PKPD_6.R")
input = read.csv("inputforMulti.csv") #load all the data
#input = read.csv("shortinput.csv") #use when testing and debugging
dose_s = as.character(unlist(input)) #get a vector of doses
dose_c = paste('c(',dose_s,')',sep='') #make all the dose to be vectors in string
dose_vec = sapply(dose_c,function(x) parse(text=x)) #evaluate the vectors in string

## combined two functions: conccurve and countfunc -- this function will make it easier to create a data table##
cure = function(dose,interval,GA){
  GAV1= -0.0114; GAmed= 28.9 #constants
  WT = unname(weight[toString(GA)]) #get weight 
  out = conccurve(dose,interval,GA) #getting concentration from all doses
  conc = out[1:nrow(out),4][-1] #nrow=1+sum(interval)/step; 4=concentration in vc (1=time (h), 2 & 3= conc. in v1 & v2)
  count = countfunc(conc) #function to get the count output
  nonzero<-count[ ,"m"]/(1+count[,"corr"]) #stochastic estimate for prob 1 cell has not gone extinct (ksr=0 approximation)
  factor<-.406*WT*4.83e8*(1+GAV1*(GA-GAmed)) #intial average number of CFU per babyvc = thvc*WT*(1+GAV1*(GA-GA_med))
  cure<-exp(-factor*nonzero) #probability every lineage has gone extinct
  return(list(cure, conc))
}

#create a dataframe. called it table
table = data.frame(Doses =character(0), E25=numeric(0), c0.1=numeric(0),#column names
                   c23.99=numeric(0),c24.1=numeric(0),c35.99=numeric(0),c36.1=numeric(0),
                   c47.99=numeric(0),c48.1=numeric(0),c71.99=numeric(0), c72.1=numeric(0),
                   E40=numeric(0), C0=numeric(0), C23.99=numeric(0),C24.1=numeric(0),C35.99=numeric(0),C36.1=numeric(0),
                   C47.99=numeric(0),C48.1=numeric(0),C71.99=numeric(0), C72.1=numeric(0),stringsAsFactors=F)

################### generating data from each set of doses ######################
for(i in 1:length(dose_vec)){#loop thru the input dosing regimens
  print(i) #just print to check R is running when we create a table
  dose = array(eval(dose_vec[i])) #convert list of dose to array
  interval=c(0,rep(12,6),4) #this is used to determine what interval correspond to the dose
  if(length(dose)==4) interval= c(0,rep(24,3),4) 
  if(length(dose)==3) interval= c(0,rep(36,2),4)
  time = cumsum(interval[-length(interval)]) 
  parameters['times'] = list(time); parameters['doses'] = list(dose); 
  cure_indices =c(dose_s[i]) #create a vector containg infor of each row. The first column is doses
  for(GA in c(25,40)){
    cure_vec = array(unlist(cure(dose,interval,GA)[1])) #probability cure
    conc = array(unlist(cure(dose,interval,GA)[2])) #concentration
    cure_indices = c(cure_indices,cure_vec[7500],conc[c(8,2399,2408,3599,3608,4799,4808,7199,7208)])
    #adding other info we want to the row
  }
  table[nrow(table)+1,] = cure_indices #add the row to the table
}
sum = as.numeric(table[,"E25"])+as.numeric(table[,"E40"]) #sum E25 and E40
sumOutput = data.frame("Doses"= table$Doses, "E25+E40"= sum)

write.csv(table, "outputMultiregimen.csv")
write.csv(sumOutput, "sumOutput.csv")

