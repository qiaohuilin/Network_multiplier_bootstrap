library(RSpectra)
library(Matrix)

setwd("~/submit/source/")
source("graphon_simulate.R")
Rcpp::sourceCpp("matmult.cpp")

setwd("~/submit/Part1_sup_error/")

set.seed(123)


n=160

# 15 sparsity levels 
tau_n_list_new=c(0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

# two graphons to choose 
input=1 #sbm 
#input=2 #smg


for(kk in 1:15){
tau_n=tau_n_list_new[kk]

if(input==1){
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
  my_smg=sparse_object(my.graphon.2, function(x) return(tau_n))
  simulate_smg=create_simulation_function(my_smg)
  
  my_smg_dense=sparse_object(my.graphon.2, function(x) return(1))
  simulate_smg_dense=create_simulation_function(my_smg_dense)
  
  A_temp=simulate_smg_dense(10000)
  degree_temp=colSums(A_temp)
  c_hat=(sum(degree_temp)/2)/(choose(10000,2))
}
  
  TT=10^6
  
  dg_hat_true=NULL
  numtri_hat_true=NULL
  twostar_hat_true=NULL
  trans_hat_true=NULL
  Sn_hat=NULL
  Sn_hat_ts=NULL
  
  for(tt in 1:TT){
    A=simulate_smg(n)
    
    #edge density
    degree_temp=colSums(A)
    dg_hat_true[tt]=sum(degree_temp/2)/(choose(n,2)*c_hat*tau_n)
    
    ### cpp function for matrix multiplication A^2
    A2=eigenMatMult(A,A)
    diag(A2)=0
    
    ##########note we can also do  ########
    ## A1=Matrix(A,sparse=T)####
    ## A2= A1 %*% A1 ####
    ## diag(A2)=0
    
    A3= A2* A   # dot product here
    num_tri_i=colSums(A3)/2
    num_tri_tot=sum(num_tri_i/3)
    numtri_hat_true[tt]=num_tri_tot/(choose(n,3)*c_hat^3*tau_n^3)
    
    
    r=3
    g1=num_tri_i/(choose(n-1,2)*c_hat^3*tau_n^3)-numtri_hat_true[tt]
    Sn_hat[tt]=sqrt(r^2/n^2*sum(g1^2))
    
    #two-stars density
    twostar_tot=sum(colSums(A2*(1-A)))/2
    twostar_hat_true[tt]= twostar_tot/(choose(n,3)*c_hat^2*tau_n^2)
    
    
    r=3
    h1ts=colSums(A2 * (1-A)) + choose(degree_temp,2)- num_tri_i
    g1ts=h1ts/(choose(n-1,2)*c_hat^2*tau_n^2)-twostar_hat_true[tt]
    Sn_hat_ts[tt]=sqrt(r^2/n^2*sum(g1ts^2))
    
    #normalized transitivity
    trans_hat_true[tt]= numtri_hat_true[tt]*3/twostar_hat_true[tt]
    print(c(input,kk,tt))
    
  }
 
  
  true=list(
    dg_hat_true=dg_hat_true,
    numtri_hat_true=numtri_hat_true,
    twostar_hat_true=twostar_hat_true,
    trans_hat_true=trans_hat_true,
    Sn_hat=Sn_hat,
    Sn_hat_ts=Sn_hat_ts)
  
  if(input=1){
    save(true,file=sprintf("truth_sbm-kk=%s-n=160.RData",kk))
  }
  if(input=2){
  save(true,file=sprintf("truth_smg-kk=%s-n=160.RData",kk))
}


}

