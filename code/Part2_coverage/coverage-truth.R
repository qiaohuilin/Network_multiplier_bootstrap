library(RSpectra)
library(Matrix)

setwd("../source/")
#load graphon simulate rcpp functions
source("graphon_simulate.R")
#load fast matrix multipication functions (optional - can use other packages)
Rcpp::sourceCpp("matmult.cpp")

setwd("../Part2_coverage/")

tau_n_list=c(0.05,0.07,0.09,0.2,0.6,1.0)

n=15000

dg_hat_true=matrix(0,2,10)
tri_hat_true=matrix(0,2,10)
twostar_hat_true=matrix(0,2,10)
trans_hat_true=matrix(0,2,10)
Sn_hat=matrix(0,2,10)
Sn_hat_ts=matrix(0,2,10)

for(input in 1:2){
  for(kk in 1:10){
  if(input ==1 ){
    tau_n= tau_n_list[kk]
    
    my_sbm_sparse=sparse_object(my.sbm, function(x) return(tau_n))
    simulate_sbm=create_simulation_function(my_sbm_sparse)
    
    my_sbm_dense=sparse_object(my.sbm, function(x) return(1))
    simulate_sbm_dense=create_simulation_function(my_sbm_dense)
    
    #estimate sparsity parameter from a large size graph
    A_temp=simulate_sbm_dense(10000)
    di=colSums(A_temp)
    c_hat=(sum(di)/2)/(choose(10000,2))
     
    A=simulate_sbm(n)
    }
    
    if(input==2){
      tau_n= tau_n_list[kk]
      my_smg=sparse_object(my.graphon.2, function(x) return(tau_n))
      simulate_smg=create_simulation_function(my_smg)
      
      my_smg_dense=sparse_object(my.graphon.2, function(x) return(1))
      simulate_smg_dense=create_simulation_function(my_smg_dense)
      
      A_temp=simulate_smg_dense(10000)
      di=colSums(A_temp)
      c_hat=(sum(di)/2)/(choose(10000,2))
      
      A=simulate_smg(n)
      }
    
    #edge density
    di=colSums(A)
    dg_hat_true[input,kk]=sum(di/2)/(choose(n,2)*c_hat*tau_n)
    
    ### cpp function for matrix multiplication A^2
    A2=eigenMatMult(A,A)
    diag(A2)=0
    
    ##########note we can also do  ########
    ## A1=Matrix(A,sparse=T)####
    ## A2= A1 %*% A1 ####
    ## diag(A2)=0
    
    A3= A2*A
    num_tri_i=colSums(A3)/2
    num_tri_tot=sum(num_tri_i/3)
    tri_hat_true[input,kk]=num_tri_tot/(choose(n,3)*c_hat^3*tau_n^3)
    
    
    r=3
    g1=num_tri_i/(choose(n-1,2)*c_hat^3*tau_n^3)-tri_hat_true[input,kk]
    Sn_hat[input,kk]=sqrt(r^2/n^2*sum(g1^2))
    
    #two-stars density
    twostar_tot=sum(colSums(A2*(1-A)))/2
    twostar_hat_true[input,kk]= twostar_tot/(choose(n,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    
    
    r=3
    h1ts=(colSums(A2 * (1-A)) + choose(di,2)- num_tri_i)/(choose(n-1,2)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    g1ts=h1ts-twostar_hat_true[input,kk]
    Sn_hat_ts[input,kk]=sqrt(r^2/n^2*sum(g1ts^2))
    
    #normalized transitivity
    trans_hat_true[input,kk]= tri_hat_true[input,kk]*3/twostar_hat_true[input,kk]
    print(c(input,kk))
    
    
    true=list(
      dg_hat_true=dg_hat_true,
      tri_hat_true=tri_hat_true,
      twostar_hat_true=twostar_hat_true,
      trans_hat_true=trans_hat_true,
      Sn_hat=Sn_hat,
      Sn_hat_ts=Sn_hat_ts)
    
    save(true,file="coverage-truth-6taus.RData")
  }
}
