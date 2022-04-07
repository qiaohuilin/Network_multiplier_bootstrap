library(RSpectra)
library(Matrix)
setwd("../source")
Rcpp::sourceCpp("matmult.cpp")

setwd("../realdata1_congress/")


data=read.csv("Hall_votes.csv")
library(dplyr)
library(tidyr)

data %>% filter(congress >80 & congress < 113) -> data
congress_list=c(81:112) #use data from 81 to 112 congress 
A_list=list(32)
n_list=vector()
threshold=c(124,76,77,80,99,93,125,117,200,257,241,340,595,712,889,690,482,518,508,563,537,516,612,763,747,786,659,781,766,1122,1094,838)
#threshold is calculated from andris2015 "The rise of partisanship and super-cooperators in the US House of Representatives"

party=read.csv("Hall_members.csv")
party= party %>% select(icpsr,party_code)
party=party %>% group_by(icpsr) %>% summarize(party_code = first(party_code))

#construct network from votes data 
#for each congress, an edge btween members are formed when the times they both vote yay exceed threshol
for(i in 1:32){
  temp=data %>% filter(congress == congress_list[i])
  tb=temp %>% select(icpsr,rollnumber,cast_code) %>% pivot_wider(names_from = rollnumber,values_from = cast_code)
  tb=left_join(tb, party,by = c("icpsr" = "icpsr"))
  
  tb=tb %>% filter(party_code==200 | party_code==100)
  #select only the two major parties 
  tb2=tb %>% select(-c("party_code","icpsr")) %>% as.matrix()
  colnames(tb2) <- NULL
  tb2[is.na(tb2)]=0
  
  tbyay=matrix(0,dim(tb2)[1],dim(tb2)[2])
  tbyay[tb2==1]=1
  tbnay=matrix(0,dim(tb2)[1],dim(tb2)[2])
  tbnay[tb2==6]=1
  
  simcomp=tbyay %*% t(tbyay) + tbnay %*% t(tbnay)
  diag(simcomp)=0
  
  A=matrix(0,dim(simcomp)[1],dim(simcomp)[2])
  A[simcomp>=threshold[i]]=1
  diag(A)=0
  A_list[[i]]=A
  n_list[i]=dim(A)[1]
}

# calculate basics for the networks constructed 
rhon_hat=vector()
degree_tot=vector()
for(i in 1:32){
  A=A_list[[i]]
  n=n_list[i]
  degree_tot[i]=sum(colSums(A))/2
  rhon_hat[i]=(sum(colSums(A))/2)/(choose(n,2))
}


############################################
#check threshold method from andris 
party=read.csv("Hall_members.csv")
party= party %>% select(icpsr,party_code)
party=party %>% group_by(icpsr) %>% summarize(party_code = first(party_code))

#check with congress 109 data -- you can check with any year's data 
temp=data %>% filter(congress == 109)
tb=temp %>% select(icpsr,rollnumber,cast_code) %>% pivot_wider(names_from = rollnumber,values_from = cast_code)
tb=left_join(tb, party,by = c("icpsr" = "icpsr"))

tb=tb %>% filter(party_code==200 | party_code==100)
tb2=tb %>% select(-c("party_code","icpsr")) %>% as.matrix()
colnames(tb2) <- NULL
tb2[is.na(tb2)]=0

tbyay=matrix(0,dim(tb2)[1],dim(tb2)[2])
tbyay[tb2==1]=1
tbnay=matrix(0,dim(tb2)[1],dim(tb2)[2])
tbnay[tb2==6]=1

simcomp=tbyay %*% t(tbyay) + tbnay %*% t(tbnay)
diag(simcomp)=0
SP1=simcomp[tb$party_code==100,tb$party_code==100]
SP2=simcomp[tb$party_code==200,tb$party_code==200]
SP=c(as.vector(SP1[SP1!=0]),as.vector(SP2[SP2!=0]))
CP=as.vector(simcomp[tb$party_code==200,tb$party_code==100])

##### use plot to check if the andris2015 threshold makes sense -- it is the point where the two histograms intersect
library(ggplot2)
ggplot() + 
  geom_density(aes(SP, col = "SP")) +
  geom_density(aes(CP, col = "CP"))

install.packages("cowplot")
library(cowplot)
theme_set(theme_minimal_grid())
ggplot() + 
  geom_histogram(aes(CP, ..density.., fill = "CP")) +
  # geom_density(aes(CP, col = "CP")) +
  NULL
ggplot() + 
  geom_histogram(aes(SP, ..density.., fill = "SP")) +
  # geom_density(aes(CP, col = "CP")) +
  NULL


#############################################
#from the networks constructed 
#calculate cross-party edge and triangles
cped <- vector()
cptri <- vector()
party_label=list()

for(i in 1:32){
  #first extract party label
  temp=data %>% filter(congress == congress_list[i])
  tb=temp %>% select(icpsr,rollnumber,cast_code) %>% pivot_wider(names_from = rollnumber,values_from = cast_code)
  tb=left_join(tb, party,by = c("icpsr" = "icpsr"))
  tb=tb %>% filter(party_code==200 | party_code==100)
  party_label[[i]]=tb$party_code
  #the network we constructed
  A=A_list[[i]]
  n=dim(A)[1]
  rhon=rhon_hat[i]
  
  #calculate cross-party edge (cped)
  cped[i]=sum(A[tb$party_code==200,tb$party_code==100])/(choose(n,2)*rhon)
  
  #calculate cross-party triangle (cptri)
  A2=eigenMatMult(A,A)
  A3= A2* 1*(A==1)
  num_tri_temp=colSums(A3)/2
  num_tri_tot=sum(num_tri_temp)/3
  Asp1=A[party_label[[i]]==100,party_label[[i]]==100]
  Asp2=A[party_label[[i]]==200,party_label[[i]]==200]
  A3sp1=eigenMatMult(Asp1,Asp1)*Asp1
  A3sp2=eigenMatMult(Asp2,Asp2)*Asp2
  numtri_sp1=colSums(A3sp1)/2
  numtri_sp2=colSums(A3sp2)/2
  numtri_cp <- vector()
  numtri_cp[party_label[[i]]==100]=num_tri_temp[party_label[[i]]==100]-numtri_sp1
  numtri_cp[party_label[[i]]==200]=num_tri_temp[party_label[[i]]==200]-numtri_sp2
  numtri_cp_tot=sum(numtri_cp)/3
  cptri[i]=numtri_cp_tot/(choose(n,3)*rhon^3)
  }

save(res=list(A_list=A_list,rhon_hat=rhon_hat,cped=cped,cptri=cptri),file='data.RData')



