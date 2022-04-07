library(RSpectra)
library(Matrix)

setwd("../source/")

source("graphon_simulate.R")
Rcpp::sourceCpp("matmult.cpp")
Rcpp::sourceCpp("brute_counts_4cyc_3pth.cpp")
Rcpp::sourceCpp("approx_counts_bts.cpp")
Rcpp::sourceCpp("approx_counts_us_original.cpp")

setwd("../Part1_sup_error/")

set.seed(123)

########### Basic parameters ###################

nlist=c(200,400,600,800,1000)
tau_n=0.1 # or any other values

input=1 #sbm 
#input=2 #SM-G

B=2000 # for bootstrap


########### Create space holders ####################
############ 5 size of graphs, 1 sparsity level, each only do 1 graph ########


tiktok_pre_linear=array(0,c(5,1))
tiktok_pre=array(0,c(5,1))
tiktok_ew=array(0,c(5,1))
tiktok_gauprod=array(0,c(5,1))
tiktok_linqua=array(0,c(5,1))
tiktok_add=array(0,c(5,1))
tiktok_lz=array(0,c(5,1))
tiktok_sub=array(0,c(5,1))
tiktok_eg=array(0,c(5,1))

fourcycle_hat=vector()

fourcycle_mb=matrix(0,1,B)


Tnhat_ut=matrix(0,5,1)
Tnhat_kato=matrix(0,5,1)

time_us=matrix(0,5,1)
time_kato=matrix(0,5,1)
tiktok_btl = matrix(0,5,1)

x=seq(-3,3,by = 0.1)
EWCDF_sd <- function(x,n,xi1,r,Eg13,Eg1g1g2){
  gx=pnorm(x)-dnorm(x)*(x^2-1)/(6*sqrt(n)*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2)
  return(gx)
}


for(kk in 1:5){
  if(input ==1 ){
    #### create sbm graphs from size n and tau_n 
    n=nlist[kk]
    tau_n=0.1
    
    my_sbm_sparse=sparse_object(my.sbm, function(x) return(tau_n))
    simulate_sbm=create_simulation_function(my_sbm_sparse)
    
    my_sbm_dense=sparse_object(my.sbm, function(x) return(1))
    simulate_sbm_dense=create_simulation_function(my_sbm_dense)
    
    #estimate sparsity parameter from a large size graph
    A_temp=simulate_sbm_dense(10000)
    degree_temp=colSums(A_temp)
    c_hat=(sum(degree_temp)/2)/(choose(10000,2))
  }
  
  if(input==2){
    #### create smg graphs from size n and tau_n
    
    n=nlist[kk]
    tau_n=0.1
    
    my_smg=sparse_object(my.graphon.2, function(x) return(tau_n))
    simulate_smg=create_simulation_function(my_smg)
    
    my_smg_dense=sparse_object(my.graphon.2, function(x) return(1))
    simulate_smg_dense=create_simulation_function(my_smg_dense)
    
    A_temp=simulate_smg_dense(10000)
    degree_temp=colSums(A_temp)
    c_hat=(sum(degree_temp)/2)/(choose(10000,2))
    }
  
  TnLboot_ut=matrix(0,1,B) #temporary for each kk
  TnLboot_kato=matrix(0,1,B)
  

  # time for one simulated graph 
  ii=1
 
  #simulate graph
  if(input==1){
    A=simulate_sbm(n)
  }
  if(input==2){
    A= simulate_smg(n)
  }
  
  #pre-calculations
  st=Sys.time()
  fourcycle_hat[ii]=count_4cycles(n,A)/(choose(n,4)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  
  r=4
  g1=count_4cycles_h1(n,A)/(choose(n-1,3)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)-fourcycle_hat[ii]
  Sn_hat_formb= sqrt(r^2/n^2*sum(g1^2))
  ed=Sys.time()
  tiktok_pre_linear[kk,ii]=as.numeric((ed-st),units="secs")
  
  g2tilde=count_4cycles_h2(n,A)/(choose(n-2,2)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)-fourcycle_hat[ii]
  diag(g2tilde)=0
  g2=g2tilde-g1*matrix(1,n,n)- t(g1*matrix(1,n,n))
  diag(g2)=0
  
  ed=Sys.time()
  tiktok_pre[kk,ii]=as.numeric((ed-st),units="secs")
  
  #Edgeworth Expansion
  st=Sys.time()
  Eg13= mean(g1^3)
  Eg1g1g2=sum(g2*(g1*matrix(1,n,n))* t(g1*matrix(1,n,n)))/(2*choose(n,2))
  xi1=sqrt(mean(g1^2))
  
  ew=EWCDF_sd(x,n,xi1,3,Eg13,Eg1g1g2)
  
  ed=Sys.time()
  tiktok_ew[kk,ii]=as.numeric((ed-st),units="secs")
  
 
  #MB-Q
  st=Sys.time()
  for(jj in 1:B){
    ee1= rnorm(n,1,sqrt(1/2)) 
    ee2= rnorm(n,1,sqrt(1/3))
    e_vec=ee1 * ee2
    
    wg2mat=(e_vec*matrix(1,n,n) * t(e_vec*matrix(1,n,n)) - matrix(1,n,n)) * h2 - (e_vec*g1*matrix(1,n,n)) - t(e_vec*g1*matrix(1,n,n)) 
    diag(wg2mat)=0
    wg2=sum(wg2mat)/2
    fourcycle_mb[ii,jj]=fourcycle_hat[ii]+r/n*sum((e_vec-1)*g1)+r*(r-1)/(n*(n-1))*wg2 
    print(c(kk,ii,jj,"MBQ"))
  }
  ed=Sys.time()
  tiktok_linqua[kk,ii]=as.numeric((ed-st),units="secs")
  
 
  #MB-L
  st=Sys.time()
  for(jj in 1:B){
    ee1= rnorm(n,1,sqrt(1/2)) 
    ee2= rnorm(n,1,sqrt(1/3))
    e_vec=ee1 * ee2
    fourcycle_mb[ii,jj]=fourcycle_hat[ii]+r/n*sum((e_vec-1)*g1)
    print(c(kk,ii,jj,"MBL"))
  }
  ed=Sys.time()
  tiktok_add[kk,ii]=as.numeric((ed-st),units="secs")
  
  
  #LS
  st=Sys.time()
  if(input ==1 ){
    eg=eigs_sym(A,20,which="LM") 
    qct= c(1:2)
    egval=eg$values[qct]
    egvec=eg$vectors[,qct]
    X=egvec %*% diag(sqrt(egval))
  }else{
    eg=eigs_sym(A,20,which="LM")
    qct=1 #eg$values > 2.1*n^(0.5)*sum(degree_temp/2)/(choose(n,2))
    egval=eg$values[qct]
    egvec=eg$vectors[,qct]
    X=egvec * sqrt(egval) 
  }
  
  P=X %*% t(X)
  diag(P)=0
  h1lz=count_4cycles_h1(n,P)/(choose(n-1,3)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  
  r=4
  fourcycle_lz=count_4cycles(n,P)/(choose(n,4)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  g1lz=h1lz-fourcycle_lz
  Sn_hat_lz= sqrt(r^2/n^2*sum(g1lz^2))
  
  W= rmultinom(B,size=n,prob = rep(1/n,n))
  for(jj in 1: B) {
    fourcycle_mb[ii,jj]=  fourcycle_lz+ 3*(sum(W[,jj]*h1lz)/n-fourcycle_lz) #3*(sum(W[,jj]*tri_i)/n) # #(sum(diag(A3))/6)/(choose(n,3)*c_hat^3*tau_n^3)
    print(c(kk,ii,jj,'LS'))
  }
  ed=Sys.time()
  tiktok_lz[kk,ii]=as.numeric((ed-st),units="secs")
  
 
  #SS
  st=Sys.time()
  Sn_hat_formbsb=NULL
  jj=1
 
  b=0.5*n
  idx=sample(1:n,b,replace = F)
  subA=A[idx,idx]
 
  fourcycle_mb[ii,jj]=count_4cycles(b,subA)/(choose(b,4)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  
  g1sb=count_4cycles_h1(b,subA)/(choose(b-1,3)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)-fourcycle_mb[ii,jj]
  Sn_hat_formbsb= sqrt(r^2/b^2*sum(g1sb^2))
  print(c(kk,ii,jj,"SS"))
  
  ed=Sys.time()
  tiktok_sub[kk,ii]=as.numeric((ed-st)*2000,units="secs")
  
  
  #EG
  st=Sys.time()
  jj=1
  idx=sample(n,n,replace = T)
  P=A[idx,idx]
  
  fourcycle_mb[ii,jj]= count_4cycles(n,P)/(choose(n,4)*c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  
  print(c(kk,ii,jj,"EG"))

  ed=Sys.time()
  tiktok_eg[kk,ii]=as.numeric((ed-st)*2000,units="secs")
  
  
  
  #######################################################
  #fast-sketched
  
  ############# using our ustatistic sampling method #########
  N=floor(50*log(n))
  r=4
  st=Sys.time()
  M=cal_4cycle_h1(n,N,A)
  h1=M[,N]/(c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)
  ed=Sys.time()
  time_us[kk,ii]=as.numeric((ed-st),units="secs")
  
  Tnhat_ut[kk,ii]=mean(h1)
  g1=h1-Tnhat_ut[kk,ii]
  Sn_hat_formb_ut= sqrt(r^2/n^2*sum(g1^2))
  
  ############## using Kato method ####################
  
  N3=floor(n*50*log(n)/r)
  
  st=Sys.time()
  h1_kato=cal_kato_4cycle_h1(n,N3,A)/(c_hat^4*tau_n^4*(1-c_hat*tau_n)^2)#-Tnhat_kato[ii,tt]
  ed=Sys.time()
  time_kato[kk,ii]=as.numeric((ed-st),units="secs")
  
  Tnhat_kato[kk,ii]=mean(h1_kato)
  g1kato=h1_kato-Tnhat_kato[kk,ii]
  Sn_hat_formb_kato= sqrt(r^2/n^2*sum(g1kato^2))
  
  
  
  #MB-L-apx
  st=Sys.time()
  for(jj in 1:B){
    ee1= rnorm(n,1,sqrt(1/2)) 
    ee2= rnorm(n,1,sqrt(1/3))
    e_vec=ee1 * ee2
    TnLboot_ut[ii,jj]=Tnhat_ut[kk,ii]+ 3/n*sum((e_vec-1)*g1)
    TnLboot_kato[ii,jj]=Tnhat_kato[kk,ii]+ 3/n*sum((e_vec-1)*g1kato)
    print(c(kk,ii,jj,"MBLappx"))
  }
  ed=Sys.time()
  tiktok_btl[kk,ii]=as.numeric((ed-st),units="secs")

  
  res=list(
    tiktok_pre_linear=tiktok_pre_linear,
    tiktok_pre=tiktok_pre,
    tiktok_ew=tiktok_ew,
    tiktok_gauprod=tiktok_gauprod,
    tiktok_linqua=tiktok_linqua,
    tiktok_add=tiktok_add,
    tiktok_lz=tiktok_lz,
    tiktok_sub=tiktok_sub,
    tiktok_eg=tiktok_eg,
    
    
    
    tiktok_btl=tiktok_btl,
    
    timeg1=timeg1,
    timeg1_kato=timeg1_kato
  
    
  )
  if(input==1){
    save(res,file="time-5n-nlist200to1000-4cycles-sbm.RData")
  }
  if(input==2){
    save(res,file="time-5n-nlist200to1000-4cycles-smg.RData")
  }
}

