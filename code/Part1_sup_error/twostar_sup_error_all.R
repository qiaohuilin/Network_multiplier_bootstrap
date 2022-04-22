library(RSpectra)
library(Matrix)

setwd("../source/")
#load graphon simulate rcpp functions
source("graphon_simulate.R")
#load fast matrix multipication functions (optional - can use other packages)
Rcpp::sourceCpp("matmult.cpp")
#load rcpp functions for approximate counts using using Chen and Kato Bernoulli tuple sampling methods
Rcpp::sourceCpp("approx_counts_bts.cpp")
#load rcpp functions for approximate counts using our original ustatistic sampling method
Rcpp::sourceCpp("approx_counts_us_original.cpp")

setwd("../Part1_sup_error/")

set.seed(123)

########### Basic parameters ###################
n=160

### 15 sparsity we used #############
tau_n_list=c(0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

input=1 #sbm 
#input=2 #SM-G

B= 2000  # for bootstrap 

########### Create space holders ####################
############ 15 sparsity level, each level 30 graphs########
############ each graph cdf (-3,3),0.1 grid, so 61 spots######### 

twostar_hat=matrix(0,15,30)

#sup errors
sup_ew=array(0,c(15,30))
sup_gauprod=array(0,c(15,30))
sup_linqua=array(0,c(15,30))
sup_add=array(0,c(15,30))
sup_lz=array(0,c(15,30))
sup_sub=array(0,c(15,30))
sup_eg=array(0,c(15,30))

#cdfs
cdfapx_true=array(0,c(15,61))
cdfapx_ew=array(0,c(15,30,61))
cdfapx_gauprod=array(0,c(15,30,61))
cdfapx_linqua=array(0,c(15,30,61))
cdfapx_add=array(0,c(15,30,61))
cdfapx_lz=array(0,c(15,30,61))
cdfapx_sub=array(0,c(15,30,61))
cdfapx_eg=array(0,c(15,30,61))

#time
tiktok_pre=array(0,c(15,30))
tiktok_ew=array(0,c(15,30))
tiktok_gauprod=array(0,c(15,30))
tiktok_linqua=array(0,c(15,30))
tiktok_add=array(0,c(15,30))
tiktok_lz=array(0,c(15,30))
tiktok_sub=array(0,c(15,30))
tiktok_eg=array(0,c(15,30))


sup_us=array(0,c(15,30))
sup_kato=array(0,c(15,30))

cdfapx_true=array(0,c(15,61))
cdfapx_us=array(0,c(15,30,61))
cdfapx_kato=array(0,c(15,30,61))

Tnhat_us=matrix(0,15,30)
Tnhat_kato=matrix(0,15,30)
time_us=matrix(0,15,30)
time_kato=matrix(0,15,30)
tiktok_btl = matrix(0,15,30)

####### function of Edgeworth expansion ###############
EWCDF_sd <- function(x,n,xi1,r,Eg13,Eg1g1g2){
  gx=pnorm(x)-dnorm(x)*(x^2-1)/(6*sqrt(n)*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2)
  return(gx)
}


for(kk in 1:15){
  if(input ==1 ){
    # this creates sbm of sparsity tau_n, kk is index of tau_n 
    
    tau_n= tau_n_list[kk]
    
    my_sbm_sparse=sparse_object(my.sbm, function(x) return(tau_n))
    simulate_sbm=create_simulation_function(my_sbm_sparse)
    
    my_sbm_dense=sparse_object(my.sbm, function(x) return(1))
    simulate_sbm_dense=create_simulation_function(my_sbm_dense)
    
    #estimate sparsity parameter from a large size graph
    A_temp=simulate_sbm_dense(10000)
    di=colSums(A_temp)
    c_hat=(sum(di)/2)/(choose(10000,2))
    
    # filename of simulated truth to load 
    filename=sprintf("simulated-truth/truth_sbm-kk=%s-n=160.RData",kk)
    
  }
  
  if(input==2){
    # this creates smg of sparsity tau_n 
    tau_n= tau_n_list[kk]
    my_smg=sparse_object(my.graphon.2, function(x) return(tau_n))
    simulate_smg=create_simulation_function(my_smg)
    
    my_smg_dense=sparse_object(my.graphon.2, function(x) return(1))
    simulate_smg_dense=create_simulation_function(my_smg_dense)
    
    A_temp=simulate_smg_dense(10000)
    di=colSums(A_temp)
    c_hat=(sum(di)/2)/(choose(10000,2))
    
    # filename of simulated truth
    filename=sprintf("simulated-truth/truth_smg-kk=%s-n=160.RData",kk)
    
  }
  
  #load simulated truth under this sparsity 
  load(filename)
  x=seq(-3,3,by = 0.1)
  trueT3=(true$twostar_hat_true-mean(true$twostar_hat_true))/mean(true$Sn_hat_ts)
  truecdf3=sapply(1:length(x),function(jj) mean(trueT3<x[jj]))
  cdfapx_true[kk,]=truecdf3
  
  #temporary spaceholder for each sparsity index kk
  # because B is big we don't want to store it for every kk. 
  twostar_mb=matrix(0,30,B)
  TnLboot_us=matrix(0,30,B) 
  TnLboot_kato=matrix(0,30,B)
  
  # each sparsity level we have 30 graphs
  for(ii in 1:30){
    
    #simulate graph
    if(input==1){
      A=simulate_sbm(n)
    }
    if(input==2){
      A= simulate_smg(n)
    }
    
    ################# #pre-calculations ######################
    st=Sys.time()
    
    ### cpp function for matrix multiplication A^2
    A2=eigenMatMult(A,A)
    diag(A2)=0
    
    ##########note we can also do  ########
    ## A1=Matrix(A,sparse=T)####
    ## A2= A1 %*% A1 ####
    ## diag(A2)=0
    
    A3= A2* A
    num_tri_i=colSums(A3)/2
    
    di=colSums(A)
  
    twostar_tot=sum(colSums(A2*(1-A)))/2
    twostar_hat[kk,ii]= twostar_tot/(choose(n,3)*c_hat^2*tau_n^2)
    
    
    r=3
    h1ts=(colSums(A2 * (1-A)) + choose(di,2)- num_tri_i)/(choose(n-1,2)*c_hat^2*tau_n^2)
    h2ts=((di*matrix(1,n,n)+t(di*matrix(1,n,n))-2-2*A3)*A+(A2*(1-A)))/(choose(n-2,1)*c_hat^2*tau_n^2)
    diag(h2ts)=0
    
    g1ts=h1ts-twostar_hat[kk,ii]
    Sn_hat_formb=sqrt(r^2/n^2*sum(g1ts^2))
    
    ######## g2(ij) is h2(ij)-Tnhat- g1(i)-g1(j) #############
    g2ts=h2ts-twostar_hat[kk,ii]-g1ts*matrix(1,n,n)- t(g1ts*matrix(1,n,n))
    diag(g2ts)=0
    
    ######## g2tilde(ij) is h2(ij)-Tnhat ###############
    g2tildets=h2ts-twostar_hat[kk,ii]
    diag(g2tildets)=0
    ed=Sys.time()
    tiktok_pre[kk,ii]=ed-st
    
    ################## Edgeworth Expansion ##################
    st=Sys.time()
    Eg13= mean(g1ts^3)
    Eg1g1g2=sum(g2ts*(g1ts*matrix(1,n,n))* t(g1ts*matrix(1,n,n)))/(2*choose(n,2))
    xi1=sqrt(mean(g1ts^2))
    
    ew=EWCDF_sd(x,n,xi1,3,Eg13,Eg1g1g2)
    
    sup_ew[kk,ii]=max(abs(ew-truecdf3))
    cdfapx_ew[kk,ii,]=ew
    
    ed=Sys.time()
    tiktok_ew[kk,ii]=ed-st
    
    
    ################## MB-M ################################
    st=Sys.time()
    for(jj in 1:B){
      
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      
      A_mb=e_vec*A
      A_mb=Matrix(A_mb,sparse = T)
      
      ex=sum(e_vec)
      ey=sum(e_vec^2)
      ez=sum(e_vec^3)
      sumeijk=(ex^3-3*(ex*ey-ez)-ez)/6
      sumeij=(ex^2-ey)/2
      
      
      WW=diag(e_vec)
      DD=WW %*% A %*% WW %*% A %*% WW
      sumtwostar=(sum(sum(DD))-sum(diag(DD))-sum(diag(DD %*% A)))/2
      sum_center_twostar=sumeijk* twostar_tot/(choose(n,3))
      
      twostar_mb[ii,jj]= twostar_tot/(choose(n,3)*c_hat^2*tau_n^2)+ (sumtwostar-sum_center_twostar)/(choose(n,3)*c_hat^2*tau_n^2)
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_gauprod[kk,ii]=ed-st
    
    Tmb_gauprod=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/Sn_hat_formb #sqrt(var(twostar_mb[ii,]))#Sn_hat_formb
    # if center at tnhat here would be
    # Tmb_gauprod=(twostar_mb[ii,]-twostar_hat[kk,ii])/Sn_hat_formb
    Tmbcdf_gauprod=sapply(1:length(x),function(qq) mean(Tmb_gauprod<x[qq]))
    sup_gauprod[kk,ii]=max(abs(Tmbcdf_gauprod-truecdf3))
    cdfapx_gauprod[kk,ii,]= Tmbcdf_gauprod
    
    ################## MB-Q ################################
    st=Sys.time()
    for(jj in 1:B){
      
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      
      wg2mat=(e_vec*matrix(1,n,n) * t(e_vec*matrix(1,n,n)) - matrix(1,n,n)) * g2tildets - (e_vec*g1ts*matrix(1,n,n)) - t(e_vec*g1ts*matrix(1,n,n)) 
      diag(wg2mat)=0
      wg2=sum(wg2mat)/2
      twostar_mb[ii,jj]=twostar_hat[kk,ii]+r/n*sum((e_vec-1)*g1ts)+r*(r-1)/(n*(n-1))*wg2 
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_linqua[kk,ii]=ed-st
    
    Tmb_linqua=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/Sn_hat_formb #sqrt(var(twostar_mb[ii,]))#Sn_hat_formb
    # if center at tnhat here would be
    # Tmb_linqua=(twostar_mb[ii,]-twostar_hat[kk,ii])/Sn_hat_formb
    Tmbcdf_linqua=sapply(1:length(x),function(qq) mean(Tmb_linqua<x[qq]))
    sup_linqua[kk,ii]=max(abs(Tmbcdf_linqua-truecdf3))
    cdfapx_linqua[kk,ii,]= Tmbcdf_linqua
    
    ################## MB-L ################################
    st=Sys.time()
    for(jj in 1:B){
      
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      
      twostar_mb[ii,jj]= twostar_hat[kk,ii]+ r/n*sum((e_vec-1)*g1ts)
      
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_add[kk,ii]=ed-st
    
    Tmb_add=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/Sn_hat_formb #sqrt(var(twostar_mb[ii,]))#Sn_hat_formb
    # if center at tnhat here would be
    # Tmb_add=(twostar_mb[ii,]-twostar_hat[kk,ii])/Sn_hat_formb
    Tmbcdf_add=sapply(1:length(x),function(qq) mean(Tmb_add<x[qq]))
    sup_add[kk,ii]=max(abs(Tmbcdf_add-truecdf3))
    cdfapx_add[kk,ii,]= Tmbcdf_add
    
    ################## LS ################################
    st=Sys.time()
    if(input ==1 ){
      eg=eigen(A)
      idx=order(abs(eg$values),decreasing = TRUE) 
      qct=sapply(1:2,function(i) which(idx==i))
      egval=eg$values[qct]
      egvec=eg$vectors[,qct]
      if(length(qct)==1){
        egvec=matrix(egvec, n, 1)
        P= (egvec * egval)  %*% t(egvec)
      }else{
        P= egvec %*% diag(egval) %*% t(egvec)
      }
    }
    if(input==2){
      eg=eigen(A)
      idx=order(abs(eg$values),decreasing = TRUE) 
      qct=abs(eg$values) > 2.01*sqrt(n*mean(A))
      egval=eg$values[qct]
      egvec=eg$vectors[,qct]
      if(sum(qct)==1){
        egvec=matrix(egvec, n, 1)
        P= (egvec * egval)  %*% t(egvec)
      }else{
        P= egvec %*% diag(egval) %*% t(egvec)
      }
    }
    
    diag(P)=0
    P[P>1]=1
    P[P<0]=0
    
    di2=colSums(P)
    P2=eigenMatMult(P,P)
    diag(P2)=0
    P3= P2* P
    num_tri_i2=colSums(P3)/2
    
    P4=P2 * (1-P)
    twostar_tot2=sum(colSums(P4))/2
    h1tslz=(colSums(P4) + choose(di2,2)- num_tri_i2)/(choose(n-1,2)*c_hat^2*tau_n^2)
    
    twostarhat2=twostar_tot2/(choose(n,3)*c_hat^2*tau_n^2)
    g1lz=h1tslz - twostarhat2
    Sn_hat_lz=sqrt(r^2/n^2*sum(g1lz^2))
    
    W= rmultinom(B,size=n,prob = rep(1/n,n))
    for(jj in 1: B) {
      twostar_mb[ii,jj]=  twostarhat2+ 3*(sum(W[,jj]*h1tslz)/n-twostarhat2) #3*(sum(W[,jj]*twostar_i)/n) #Untwo+ 3*(sum(W[,jj]*twostar_i)/n-Untwo) #twostar_tot/(choose(n,3)*c_hat^2*tau_n^2)
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_lz[kk,ii]=ed-st
    
    Tmb_lz=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/Sn_hat_lz #sqrt(var(twostar_mb[ii,]))#Sn_hat_formb
    # if center at tnhat here would be
    # Tmb_lz=(twostar_mb[ii,]-twostar_hat[kk,ii])/Sn_hat_lz 
    Tmbcdf_lz=sapply(1:length(x),function(qq) mean(Tmb_lz<x[qq]))
    sup_lz[kk,ii]=max(abs(Tmbcdf_lz-truecdf3))
    cdfapx_lz[kk,ii,]= Tmbcdf_lz
    
    ################## SS ################################
    st=Sys.time()
    Sn_hat_formbsb=NULL
    for(jj in 1:B){
      b=0.5*n
      idx=sample(1:n,b,replace = F)
      subA=A[idx,idx]
      subA2=eigenMatMult(subA,subA)
      diag(subA2)=0
      disb=colSums(subA)
      
      subA3= subA2* subA
      num_tri_isb=colSums(subA3)/2
      
      twostar_totsb=sum(colSums(subA2*(1-subA)))/2
      twostar_mb[ii,jj]=twostar_totsb/(choose(b,3)*c_hat^2*tau_n^2)
      
      h1tssb=(colSums(subA2 * (1-subA)) + choose(disb,2)- num_tri_isb)/(choose(b-1,2)*c_hat^2*tau_n^2)
      g1tssb=h1tssb-twostar_mb[ii,jj]
      Sn_hat_formbsb[jj]=sqrt(r^2/b^2*sum(g1tssb^2))
      
      print(c(input,kk,ii))
    }
    ed=Sys.time()
    tiktok_sub[kk,ii]=ed-st
    
    Tmb_sub=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/(Sn_hat_formbsb*sqrt(b/n))  
    # if center at tnhat here would be
    # Tmb_sub=(twostar_mb[ii,]-twostar_hat[kk,ii])/(Sn_hat_formbsb*sqrt(b/n))  
    Tmbcdf_sub=sapply(1:length(x),function(qq) mean(Tmb_sub<x[qq]))
    sup_sub[kk,ii]=max(abs(Tmbcdf_sub-truecdf3))
    cdfapx_sub[kk,ii,]= Tmbcdf_sub
    
    ################## EG ################################
    st=Sys.time()
    for(jj in 1:B){
      idx=sample(n,n,replace = T)
      P=A[idx,idx]
      di2=colSums(P)
      P2=eigenMatMult(P,P)
      diag(P2)=0
      
      twostar_tot2=sum(colSums(P2*(1-P)))/2
      twostar_mb[ii,jj]= twostar_tot2/(choose(n,3)*c_hat^2*tau_n^2)
      
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_eg[kk,ii]=ed-st
    
    Tmb_eg=(twostar_mb[ii,]-mean(twostar_mb[ii,]))/Sn_hat_formb 
    # if center at tnhat here would be
    #Tmb_eg=(twostar_mb[ii,]-twostar_hat[kk,ii])/Sn_hat_formb 
    Tmbcdf_eg=sapply(1:length(x),function(qq) mean(Tmb_eg<x[qq]))
    sup_eg[kk,ii]=max(abs(Tmbcdf_eg-truecdf3))
    cdfapx_eg[kk,ii,]= Tmbcdf_eg
    
    
    #######################################################
    #fast-sketched
    
    ############# using our ustatistic sampling method #########
    N=floor(50*log(n))
    st=Sys.time()
    M=cal_vstar_h1(n,N,A)  #see source code function
    h1=M[,N]/(c_hat^2*tau_n^2)
    ed=Sys.time()
    time_us[kk,ii]=ed-st
    
    Tnhat_us[kk,ii]=mean(h1)
    g1=h1-Tnhat_us[kk,ii]
    Sn_hat_formb_us= sqrt(r^2/n^2*sum(g1^2))
    
    ############## using Kato method ####################
    N3=floor(n*50*log(n)/r)
    
    st=Sys.time()
    h1_kato=cal_kato_vstar_h1(n,N3,A)/(choose(n-1,2)*c_hat^2*tau_n^2)
    #see source code function
    ed=Sys.time()
    time_kato[kk,ii]=ed-st
    
    Tnhat_kato[kk,ii]=mean(h1_kato)
    g1kato=h1_kato-Tnhat_kato[kk,ii]
    Sn_hat_formb_kato= sqrt(r^2/n^2*sum(g1kato^2))
    
    
    #MB-L-apx
    st=Sys.time()
    for(jj in 1:B){
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      TnLboot_us[ii,jj]=Tnhat_us[kk,ii]+ 3/n*sum((e_vec-1)*g1)
      TnLboot_kato[ii,jj]=Tnhat_kato[kk,ii]+ 3/n*sum((e_vec-1)*g1kato)
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_btl[kk,ii]=ed-st
    
    Tmb_us=(TnLboot_us[ii,]-mean(TnLboot_us[ii,]))/Sn_hat_formb_us 
    # if center at tnhat here would be
    # Tmb_us=(TnLboot_us[ii,]-Tnhat_us[kk,ii])/Sn_hat_formb_us 
    Tmbcdf_us=sapply(1:length(x),function(qq) mean(Tmb_us<x[qq]))
    sup_us[kk,ii]=max(abs(Tmbcdf_us-truecdf3))
    cdfapx_us[kk,ii,]= Tmbcdf_us
    
    Tmb_kato=(TnLboot_kato[ii,]-mean(TnLboot_kato[ii,]))/Sn_hat_formb_kato 
    # if center at tnhat here would be
    #Tmb_kato=(TnLboot_kato[ii,]-Tnhat_kato[kk,ii])/Sn_hat_formb_kato 
    Tmbcdf_kato=sapply(1:length(x),function(qq) mean(Tmb_kato<x[qq]))
    sup_kato[kk,ii]=max(abs(Tmbcdf_kato-truecdf3))
    cdfapx_kato[kk,ii,]= Tmbcdf_kato
    
    res=list(
      
      ###### sups are the ones we care about
      sup_ew=sup_ew,
      sup_gauprod=sup_gauprod,
      sup_linqua=sup_linqua,
      sup_add=sup_add,
      sup_lz=sup_lz,
      sup_sub=sup_sub,
      sup_eg=sup_eg,
      sup_us=sup_us,
      sup_kato=sup_kato,
      
      # the rest is just in case
      cdfapx_true=cdfapx_true,
      cdfapx_ew=cdfapx_ew,
      cdfapx_gauprod=cdfapx_gauprod,
      cdfapx_linqua=cdfapx_linqua,
      cdfapx_add=cdfapx_add,
      cdfapx_lz=cdfapx_lz,
      cdfapx_sub=cdfapx_sub,
      cdfapx_eg=cdfapx_eg,
      cdfapx_us=cdfapx_us,
      cdfapx_kato=cdfapx_kato,
      
      tiktok_pre=tiktok_pre,
      tiktok_ew=tiktok_ew,
      tiktok_gauprod=tiktok_gauprod,
      tiktok_linqua=tiktok_linqua,
      tiktok_add=tiktok_add,
      tiktok_lz=tiktok_lz,
      tiktok_sub=tiktok_sub,
      tiktok_eg=tiktok_eg,
      tiktok_btl=tiktok_btl,
      tiktok_btl=tiktok_btl,
      
      time_us=time_us,
      time_kato=time_kato,
      Tnhat_us=Tnhat_us,
      Tnhat_kato=Tnhat_kato
      
    )
    if(input==1){
    save(res,file="EW-MB-twostar-tau-sbm-15taus.RData")
    }
    if(input==2){
    save(res,file="EW-MB-twostar-tau-smg-15taus.RData")
    }
  }
  
}
