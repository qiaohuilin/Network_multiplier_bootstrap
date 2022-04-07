library(RSpectra)
library(Matrix)

setwd("../source/")
#load graphon simulate rcpp functions
source("graphon_simulate.R")
#load fast matrix multipication functions (optional - can use other packages)
Rcpp::sourceCpp("matmult.cpp")
#load rcpp functions for approximate counts using Chen and Kato Bernoulli tuple sampling methods
Rcpp::sourceCpp("approx_counts_bts.cpp")
#load rcpp functions for approximate counts using our original ustatistic sampling method
Rcpp::sourceCpp("approx_counts_us_original.cpp")

setwd("../Part2_coverage/")

set.seed(123)

########### Basic parameters ###################
n=500

tau_n_list=c(0.05,0.07,0.09,0.2,0.6,1.0)

input=1 #sbm 
#input=2 #SM-G

B=1000   #for bootstrap


########### Create space holders ####################
############ 6 sparsity level, each level 200 graphs########

tiktok_pre=array(0,c(6,200))
tiktok_ew=array(0,c(6,200))
tiktok_gauprod=array(0,c(6,200))
tiktok_linqua=array(0,c(6,200))
tiktok_add=array(0,c(6,200))
tiktok_lz=array(0,c(6,200))
tiktok_sub=array(0,c(6,200))
tiktok_eg=array(0,c(6,200))


dg_hat=matrix(0,6,200)
tri_hat=matrix(0,6,200)
twostar_hat=matrix(0,6,200)
trans_hat=matrix(0,6,200)


Tnhat_us=matrix(0,6,200)
Tnhat_kato=matrix(0,6,200)
time_us=matrix(0,6,200)
time_kato=matrix(0,6,200)
tiktok_btl = matrix(0,6,200)

Tnhat_us_ts=matrix(0,6,200)
Tnhat_kato_ts=matrix(0,6,200)

Tnhat_us_trans=matrix(0,6,200)
Tnhat_kato_trans=matrix(0,6,200)

var_jk_tilda_trans=matrix(0,6,200)

##### place for bootstrap samples ###
tri_bt_ew=array(0,dim=c(6,200,B))
tri_bt_linqua=array(0,dim=c(6,200,B))
tri_bt_add=array(0,dim=c(6,200,B))
tri_bt_lz=array(0,dim=c(6,200,B))
tri_bt_sub=array(0,dim=c(6,200,B))
tri_bt_eg=array(0,dim=c(6,200,B))
tri_bt_us=array(0,dim=c(6,200,B))
tri_bt_kato=array(0,dim=c(6,200,B))

ts_bt_ew=array(0,dim=c(6,200,B))
ts_bt_linqua=array(0,dim=c(6,200,B))
ts_bt_add=array(0,dim=c(6,200,B))
ts_bt_lz=array(0,dim=c(6,200,B))
ts_bt_sub=array(0,dim=c(6,200,B))
ts_bt_eg=array(0,dim=c(6,200,B))
ts_bt_us=array(0,dim=c(6,200,B))
ts_bt_kato=array(0,dim=c(6,200,B))

trans_bt_ew=array(0,dim=c(6,200,B))
trans_bt_linqua=array(0,dim=c(6,200,B))
trans_bt_add=array(0,dim=c(6,200,B))
trans_bt_lz=array(0,dim=c(6,200,B))
trans_bt_sub=array(0,dim=c(6,200,B))
trans_bt_eg=array(0,dim=c(6,200,B))
trans_bt_us=array(0,dim=c(6,200,B))
trans_bt_kato=array(0,dim=c(6,200,B))


### these are quantities for corrections. crtmat is the final correction amount ######
tri_p1mat=matrix(0,6,200)
tri_q1mat=matrix(0,6,200)
tri_crtmat=matrix(0,6,200)
tri_xi1mat=matrix(0,6,200)
tri_Eg13mat=matrix(0,6,200)
tri_Eg1g1g2mat=matrix(0,6,200)

ts_p1mat=matrix(0,6,200)
ts_q1mat=matrix(0,6,200)
ts_crtmat=matrix(0,6,200)
ts_xi1mat=matrix(0,6,200)
ts_Eg13mat=matrix(0,6,200)
ts_Eg1g1g2mat=matrix(0,6,200)

#p1 and q1 are symmetric re x, z_alpha and z_1-alpha have same p1 q1 value
p1 <- function(x,xi1,r,Eg13,Eg1g1g2){
  return(-(x^2-1)/(6*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2))
}


q1 <- function(x,xi1,r,Eg13,Eg1g1g2){
  return(1/(xi1^3)*((2*x^2+1)/6*Eg13+(r-1)*(x^2+1)/2*Eg1g1g2))
}

p1_trans <- function(x,sigma,A1,A2){
  return(-(A1*sigma^(-1)+1/6*A2*sigma^(-3)*(x^2-1)))
}

q1_trans <- function(x,B1,B2){
  return(-(B1+1/6*B2*(x^2-1)))
}

p1trans_mat=matrix(0,6,200)
q1trans_mat=matrix(0,6,200)
trans_crtmat=matrix(0,6,200)

EWCDF_sd <- function(x,n,xi1,r,Eg13,Eg1g1g2){
  gx=pnorm(x)-dnorm(x)*(x^2-1)/(6*sqrt(n)*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2)
  return(gx)
}

# EW for standardize and studentize transitivity
EWCDF_sd_trans <- function(x,n,sigma,A1,A2){
  gx=pnorm(x)+n^(-1/2)*dnorm(x)*p1_trans(x,sigma,A1,A2)
  return(gx)
}

EWCDF_stu_trans <- function(x,n,B1,B2){
  gx=pnorm(x)+n^(-1/2)*dnorm(x)*q1_trans(x,B1,B2)
  return(gx)
}



for(kk in 1:6){
  if(input ==1 ){
    tau_n= tau_n_list[kk]
    
    my_sbm_sparse=sparse_object(my.sbm, function(x) return(tau_n))
    simulate_sbm=create_simulation_function(my_sbm_sparse)
    
    my_sbm_dense=sparse_object(my.sbm, function(x) return(1))
    simulate_sbm_dense=create_simulation_function(my_sbm_dense)
    
    #estimate sparsity parameter from a large size graph
    A_temp=simulate_sbm_dense(20000)
    di=colSums(A_temp)
    c_hat=(sum(di)/2)/(choose(20000,2))
   
  }
  
  if(input==2){
    tau_n= tau_n_list[kk]
    my_smg=sparse_object(my.graphon.2, function(x) return(tau_n))
    simulate_smg=create_simulation_function(my_smg)
    
    my_smg_dense=sparse_object(my.graphon.2, function(x) return(1))
    simulate_smg_dense=create_simulation_function(my_smg_dense)
    
    A_temp=simulate_smg_dense(20000)
    di=colSums(A_temp)
    c_hat=(sum(di)/2)/(choose(20000,2))
  
  }
  
  #temporary for each kk to save space 
  dg_mb=matrix(0,200,B)
  tri_mb=matrix(0,200,B)
  twostar_mb=matrix(0,200,B)
  trans_mb=matrix(0,200,B)
  
  TnLboot_us=matrix(0,200,B)
  TnLboot_kato=matrix(0,200,B)
  TnLboot_us_ts=matrix(0,200,B) 
  TnLboot_kato_ts=matrix(0,200,B)
  TnLboot_us_trans=matrix(0,200,B) 
  TnLboot_kato_trans=matrix(0,200,B)
  
  for(ii in 1:200){
    #simulate graph
    if(input==1){
      A=simulate_sbm(n)
    }
    if(input==2){
      A= simulate_smg(n)
    }
    
    ################# #pre-calculations ######################
    st=Sys.time()
    di=colSums(A)
    A2=eigenMatMult(A,A)
    diag(A2)=0
    
    ##########note we can also do  ########
    ## A1=Matrix(A,sparse=T)####
    ## A2= A1 %*% A1 ####
    ## diag(A2)=0
    
    ########  for triangles  ##########################
    A3= A2*A
    num_tri_i=colSums(A3)/2
    tri_hat[kk,ii]= (sum(num_tri_i)/3)/(choose(n,3)*c_hat^3*tau_n^3)
    
    r=3
    h1=num_tri_i/(choose(n-1,2)*c_hat^3*tau_n^3)
    g1=h1-tri_hat[kk,ii]
    Sn_hat_formb= sqrt(r^2/n^2*sum(g1^2))
    ######## g2(ij) is h2(ij)-Tnhat- g1(i)-g1(j) #############
    h2=A3/(choose(n-2,1)*c_hat^3*tau_n^3)
    g2=h2-tri_hat[kk,ii]-g1*matrix(1,n,n)- t(g1*matrix(1,n,n))
    diag(g2)=0
    ######## g2tilde(ij) is h2(ij)-Tnhat ###############
    g2tilde=h2-tri_hat[kk,ii]
    diag(g2tilde)=0
    
    ######## for two-stars #######
    twostar_tot=sum(colSums(A2*(1-A)))/2
    twostar_hat[kk,ii]= twostar_tot/(choose(n,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    
    num_ts_i=colSums(A2 * (1-A)) + choose(di,2)- num_tri_i
    h1ts=num_ts_i/(choose(n-1,2)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    h2ts=((di*matrix(1,n,n)+t(di*matrix(1,n,n))-2-2*A3)*A+(A2*(1-A)))/(choose(n-2,1)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    diag(h2ts)=0
    
    g1ts=h1ts-twostar_hat[kk,ii]
    Sn_hat_formb_ts=sqrt(r^2/n^2*sum(g1ts^2))
    g2ts=h2ts-twostar_hat[kk,ii]-g1ts*matrix(1,n,n)- t(g1ts*matrix(1,n,n))
    diag(g2ts)=0
    g2tildets=h2ts-twostar_hat[kk,ii]
    diag(g2tildets)=0
    
    ### transitivity ##############
    trans_hat[kk,ii]=tri_hat[kk,ii]*3/twostar_hat[kk,ii]
    
    #do jackknife on the whold graph
    degree_jk=((sum(di)-2*di)/2)/(choose(n-1,2)*c_hat*tau_n)
    num_tri_jk= ((sum(num_tri_i)-3*num_tri_i)/3)/(choose(n-1,3)*c_hat^3*tau_n^3)
    twostar_jk= (twostar_tot-num_ts_i)/(choose(n-1,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    trans_jk= num_tri_jk*3/twostar_jk
    
    #calculate varjack
    var_jk_tilda_trans[kk,ii]=sum((trans_jk-mean(trans_jk))^2)
    Sn_trans_jk=sqrt(var_jk_tilda_trans[kk,ii])
    
    ed=Sys.time()
    tiktok_pre[kk,ii]=ed-st
  
    
    #Edgeworth Expansion for tri
    st=Sys.time()
    Eg13= mean(g1^3)
    Eg1g1g2=sum(g2*(g1*matrix(1,n,n))* t(g1*matrix(1,n,n)))/(2*choose(n,2))
    xi1=sqrt(mean(g1^2))
    
    
    tri_p1mat[kk,ii]=p1(1.96,xi1,3,Eg13,Eg1g1g2)
    tri_q1mat[kk,ii]=q1(1.96,xi1,3,Eg13,Eg1g1g2)
    tri_crtmat[kk,ii]=n^(-1)*Sn_hat_formb*( tri_p1mat[kk,ii]+ tri_q1mat[kk,ii])
    tri_Eg13mat[kk,ii]=Eg13
    tri_Eg1g1g2mat[kk,ii]=Eg1g1g2
    tri_xi1mat[kk,ii]=xi1
    
    #Edgeworth Expansion for two-stars
    #st=Sys.time()
    Eg13= mean(g1ts^3)
    Eg1g1g2=sum(g2ts*(g1ts*matrix(1,n,n))* t(g1ts*matrix(1,n,n)))/(2*choose(n,2))
    xi1=sqrt(mean(g1ts^2))
    
 
    ts_p1mat[kk,ii]=p1(1.96,xi1,3,Eg13,Eg1g1g2)
    ts_q1mat[kk,ii]=q1(1.96,xi1,3,Eg13,Eg1g1g2)
    ts_crtmat[kk,ii]=n^(-1)*Sn_hat_formb_ts*( ts_p1mat[kk,ii]+ ts_q1mat[kk,ii])
    ts_Eg13mat[kk,ii]=Eg13
    ts_Eg1g1g2mat[kk,ii]=Eg1g1g2
    ts_xi1mat[kk,ii]=xi1
    
    ed=Sys.time()
    tiktok_ew[kk,ii]=ed-st
    
    ########################## correction calculation ############
    a1mat=rep(0,5)
    a1mat[1]=3/twostar_hat[kk,ii]
    a1mat[2]=-3*tri_hat[kk,ii]/(twostar_hat[kk,ii]^2)
    a1mat[3]=0
    a1mat[4]=0
    a1mat[5]=0
    
    mu12mat=matrix(0,2,2)
    mu12mat[1,2]=1/n*sum(r^2*g1*g1ts)
    mu12mat[2,1]=1/n*sum(r^2*g1*g1ts)
    mu12mat[1,1]=1/n*sum(r^2*g1*g1)
    mu12mat[2,2]=1/n*sum(r^2*g1ts*g1ts)
    
    a12mat=matrix(0,5,5)
    a12mat[1,2]=-3/(twostar_hat[kk,ii]^2)
    a12mat[2,1]=-3/(twostar_hat[kk,ii]^2)
    a12mat[1,1]=0
    a12mat[2,2]=6*tri_hat[kk,ii]/(twostar_hat[kk,ii]^3)
    a12mat[3:5,3:5]=0
    
    sigma2hat=0
    A1hat=0
    for(i in 1:2){
      for(j in 1:2){
        sigma2hat=sigma2hat+ a1mat[i]* a1mat[j]*mu12mat[i,j]
        A1hat=A1hat + 1/2 * a12mat[i,j]*mu12mat[i,j]
      }
    }
    
    sigmahat=sqrt(sigma2hat)
    
    s2=0
    for(i in 1:2){
      for(j in 1:2){
        for(k in 1:2){
          for(t in 1:2){
            s2=s2+3*a1mat[i]*a1mat[j]*a12mat[k,t]*mu12mat[i,k]*mu12mat[j,t]
          }
        }
      }
    }
    s1=0
    l=rbind(3*g1*1,3*g1ts) 
    for(i in 1:2){
      for(j in 1:2){
        for(k in 1:2){
          s1=s1+a1mat[i]*a1mat[j]*a1mat[k]* 1/n*sum(l[i,]*l[j,]*l[k,])
        }
      }
    }
    s3=0
    lg1=rbind(g1,g1ts)
    lg2=list()
    lg2[[1]]=g2
    lg2[[2]]=g2ts
    rs=c(3,3)
    for(i in 1:2){
      for(j in 1:2){
        for(k in 1:2){
          s3=s3+3*(rs[i]*rs[j]*rs[k]*(rs[k]-1))*a1mat[i]*a1mat[j]*a1mat[k]*sum(lg2[[k]]*(lg1[i,]*matrix(1,n,n))* t(lg1[j,]*matrix(1,n,n)))/(2*choose(n,2))
        }
      }
    }
    A2hat=s1+s2+s3
    
  
    
    Euab=matrix(0,2,5)
    Euab[1,1]=r^2*mean(g1*g1)/n
    Euab[1,2]=r^2*mean(g1*g1ts)/n
    Euab[2,1]=r^2*mean(g1*g1ts)/n
    Euab[2,2]=r^2*mean(g1ts*g1ts)/n
    
    Euab[1,4]=r^3/n*mean(g1^3) + 2*r^3*(r-1)/n*sum(g1*matrix(1,n,n) * g2 * t(g1*matrix(1,n,n)))/(n*(n-1)) +
      2*r^2*tri_hat[kk,ii]*r^2/n*mean(g1^2)
    Euab[2,5]=r^3/n*mean(g1ts^3) + 2*r^3*(r-1)/n*sum(g1ts*matrix(1,n,n) * g2ts * t(g1ts*matrix(1,n,n)))/(n*(n-1)) + 
      2*r^2*twostar_hat[kk,ii]*r^2/n*mean(g1ts^2)
    
    Euab[1,5]=r^3/n*mean(g1*g1ts^2) + 2*r^3*(r-1)/n*sum(g1*matrix(1,n,n) * g2ts * t(g1ts*matrix(1,n,n)))/(n*(n-1)) +
      2*r^2*twostar_hat[kk,ii]*r^2/n*mean(g1*g1ts)
    Euab[2,4]=r^3/n*mean(g1ts*g1^2)+2*r^3*(r-1)/n*sum(g1ts*matrix(1,n,n) * g2 * t(g1*matrix(1,n,n)))/(n*(n-1)) +
      2*r^2*tri_hat[kk,ii]*r^2/n*mean(g1*g1ts)
    
    Euab[1,3]=r^3/n*mean(g1^2*g1ts) + r^3*(r-1)/n*sum(g1*matrix(1,n,n) * g2ts * t(g1*matrix(1,n,n)))/(n*(n-1)) +  
      r^3*(r-1)/n*sum(g1*matrix(1,n,n) * g2 * t(g1ts*matrix(1,n,n)))/(n*(n-1)) + 
      r^2*tri_hat[kk,ii]*r^2/n*mean(g1*g1ts) + 
      r^2*twostar_hat[kk,ii] * r^2/n * mean(g1^2)
    
    Euab[2,3]=r^3/n*mean(g1*g1ts^2) + r^3*(r-1)/n*sum(g1ts*matrix(1,n,n) * g2ts * t(g1*matrix(1,n,n)))/(n*(n-1)) + 
      r^3*(r-1)/n*sum(g1ts*matrix(1,n,n) * g2 * t(g1ts*matrix(1,n,n)))/(n*(n-1)) +
      r^2*twostar_hat[kk,ii]*r^2/n*mean(g1*g1ts) +
      r^2*tri_hat[kk,ii]*r^2/n*mean(g1ts^2)
    
    
    muvec=c(tri_hat[kk,ii],twostar_hat[kk,ii],r^2*mean(h1*h1ts),
            r^2*mean(h1^2),r^2*mean(h1ts^2))
    muab=matrix(0,2,2)
    muab[1,1]=muvec[4]-r^2*muvec[1]^2
    muab[1,2]=muvec[3]-r^2*muvec[1]*muvec[2]
    muab[2,1]=muvec[3]-r^2*muvec[1]*muvec[2]
    muab[2,2]=muvec[5]-r^2*muvec[2]^2
    
    cmat=rep(0,5)
    cmat[1]=2*(a12mat[1,1]*a1mat[1]*muab[1,1] + a12mat[1,1]*a1mat[2]*muab[1,2] 
               + a12mat[2,1]*a1mat[1]*muab[2,1] +  a12mat[2,1]*a1mat[2]*muab[2,2]) - 2*(
                 a1mat[1]^2)*tri_hat[kk,ii]*r^2 - 2*a1mat[1]*a1mat[2]*twostar_hat[kk,ii]*r^2
    
    cmat[2]=2*(a12mat[1,2]*a1mat[1]*muab[1,1] + a12mat[1,2]*a1mat[2]*muab[1,2] 
               + a12mat[2,2]*a1mat[1]*muab[2,1] +  a12mat[2,2]*a1mat[2]*muab[2,2]) - 2*(
                 a1mat[2]*a1mat[1])*tri_hat[kk,ii]*r^2 - 2*a1mat[2]*a1mat[2]*twostar_hat[kk,ii]*r^2
    
    cmat[3]=2*a1mat[1]*a1mat[2]
    cmat[4]=a1mat[1]*a1mat[1]
    cmat[5]=a1mat[2]*a1mat[2]
   
    s2=0
    for(a in 1:2){
      for(b in 1:5){
        s2=s2+a1mat[a]*cmat[b]*Euab[a,b]
      }
    }
    
    B1hat = A1hat * sigmahat^(-1) - 1/2*sigmahat^(-3)*n* s2
    B2hat= 6*(B1hat - A1hat * sigmahat^(-1) + 1/6* A2hat * sigmahat^(-3))
    
    
    p1trans_mat[kk,ii]=p1_trans(1.96,sigmahat,A1hat,A2hat)
    q1trans_mat[kk,ii]=q1_trans(1.96,B1hat,B2hat)
    trans_crtmat[kk,ii]=n^(-1)*sigmahat^(-1)*( p1trans_mat[kk,ii]+ q1trans_mat[kk,ii])
    
    ###########################################################
    
    ################### MB-Q ##################################
    st=Sys.time()
    for(jj in 1:B){
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      
      wg2mat=(e_vec*matrix(1,n,n) * t(e_vec*matrix(1,n,n)) - matrix(1,n,n)) * g2tilde - (e_vec*g1*matrix(1,n,n)) - t(e_vec*g1*matrix(1,n,n)) 
      diag(wg2mat)=0
      wg2=sum(wg2mat)/2
      tri_mb[ii,jj]=tri_hat[kk,ii]+r/n*sum((e_vec-1)*g1)+r*(r-1)/(n*(n-1))*wg2 
      
      
      wg2mat=(e_vec*matrix(1,n,n) * t(e_vec*matrix(1,n,n)) - matrix(1,n,n)) * g2tildets - (e_vec*g1ts*matrix(1,n,n)) - t(e_vec*g1ts*matrix(1,n,n)) 
      diag(wg2mat)=0
      wg2=sum(wg2mat)/2
      twostar_mb[ii,jj]=twostar_hat[kk,ii]+r/n*sum((e_vec-1)*g1ts)+r*(r-1)/(n*(n-1))*wg2 
      
      trans_mb[ii,jj]=tri_mb[ii,jj]*3/ twostar_mb[ii,jj]
      print(c(kk,ii,jj,"MBQ"))
    }
    ed=Sys.time()
    tiktok_linqua[kk,ii]=ed-st
    
  
    tri_bt_linqua[kk,ii,]=tri_mb[ii,]
    ts_bt_linqua[kk,ii,]=twostar_mb[ii,]
    trans_bt_linqua[kk,ii,]=trans_mb[ii,]
    
    
    
    ################### MB-L ##################################
    st=Sys.time()
    for(jj in 1:B){
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      
      
      tri_mb[ii,jj]= tri_hat[kk,ii]+r/n*sum((e_vec-1)*g1)
      twostar_mb[ii,jj]=  twostar_hat[kk,ii]+ r/n*sum((e_vec-1)*g1ts)
      trans_mb[ii,jj]= tri_mb[ii,jj]*3/twostar_mb[ii,jj]
      print(c(kk,ii,jj,"MBL"))
    }
    ed=Sys.time()
    tiktok_add[kk,ii]=ed-st
    

    tri_bt_add[kk,ii,]=tri_mb[ii,]
    ts_bt_add[kk,ii,]=twostar_mb[ii,]
    trans_bt_add[kk,ii,]=trans_mb[ii,]
    
    
    ################### LS ##################################
    st=Sys.time()
    if(input ==1 ){
      eg=eigs_sym(A,20,which="LM")
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
      eg=eigs_sym(A,20,which="LM")
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
    
    #### LS for triangles 
    P2=eigenMatMult(P,P)
    diag(P2)=0
    P3= P2* P
    num_tri_i2=colSums(P3)/2
    h1_lz=num_tri_i2/(choose(n-1,2)*c_hat^3*tau_n^3)
    trihat_lz=(sum(num_tri_i2)/3)/(choose(n,3)*c_hat^3*tau_n^3)
    
    r=3
    g1lz=h1_lz-trihat_lz
    Sn_hat_lz= sqrt(r^2/n^2*sum(g1lz^2))
    
    ### LS for twostar
    P4=P2 * (1-P)
    twostar_tot2=sum(colSums(P4))/2
    num_ts_i2=colSums(P4) + choose(di2,2)- num_tri_i2
    h1tslz=num_ts_i2 /(choose(n-1,2)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    tshat_lz=twostar_tot2/(choose(n,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    transhat2=trihat_lz*3/tshat_lz
    
    g1tslz=h1tslz - tshat_lz
    Sn_hat_tslz=sqrt(r^2/n^2*sum(g1tslz^2))
    
    #### JK var for LS ###
    num_tri_jklz= ((sum(num_tri_i2)-3*num_tri_i2)/3)/(choose(n-1,3)*c_hat^3*tau_n^3)
    twostar_jklz= (twostar_tot2-num_ts_i2)/(choose(n-1,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    trans_jklz= num_tri_jklz*3/twostar_jklz
    Sn_trans_jklz= sqrt(sum((trans_jklz-mean(trans_jklz))^2))
    
    
    W= rmultinom(B,size=n,prob = rep(1/n,n))
    for(jj in 1: B) {
      tri_mb[ii,jj]=  trihat_lz+ 3*(sum(W[,jj]*h1_lz)/n-tri_lz) #3*(sum(W[,jj]*tri_i)/n) # #(sum(diag(A3))/6)/(choose(n,3)*c_hat^3*tau_n^3)
      twostar_mb[ii,jj]=  tshat_lz+ 3*(sum(W[,jj]*h1tslz)/n-tshat_lz) #3*(sum(W[,jj]*twostar_i)/n) #Untwo+ 3*(sum(W[,jj]*twostar_i)/n-Untwo) #twostar_tot/(choose(n,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
      trans_mb[ii,jj]=3*tri_mb[ii,jj] /  twostar_mb[ii,jj] 
      print(c(kk,ii,jj,"LS"))
    }
    ed=Sys.time()
    tiktok_lz[kk,ii]=ed-st

    tri_bt_lz[kk,ii,]=tri_mb[ii,]
    ts_bt_lz[kk,ii,]=twostar_mb[ii,]
    trans_bt_lz[kk,ii,]=trans_mb[ii,]
    
    ################### SS ##################################
    st=Sys.time()
    Sn_hat_formbsb=NULL
    Sn_hat_formb_tssb=NULL
    for(jj in 1:B){
      b=0.5*n
      idx=sample(1:n,b,replace = F)
      subA=A[idx,idx]
      subA2=eigenMatMult(subA,subA)
      diag(subA2)=0
      disb=colSums(subA)
      
      subA3= subA2* subA
      num_tri_isb=colSums(subA3)/2
      tri_mb[ii,jj]=(sum(num_tri_isb)/3)/(choose(b,3)*c_hat^3*tau_n^3)
      
      g1sb=num_tri_isb/(choose(b-1,2)*c_hat^3*tau_n^3)-tri_mb[ii,jj]
      Sn_hat_formbsb= sqrt(r^2/b^2*sum(g1sb^2))
      
      
      twostar_totsb=sum(colSums(subA2*(1-subA)))/2
      twostar_mb[ii,jj]=twostar_totsb/(choose(b,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
      
      r=3
      h1tssb=(colSums(subA2 * (1-subA)) + choose(disb,2)- num_tri_isb)/(choose(b-1,2)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
      g1tssb=h1tssb-twostar_mb[ii,jj]
      Sn_hat_formb_tssb[jj]=sqrt(r^2/b^2*sum(g1tssb^2))
      
      trans_mb[ii,jj]= tri_mb[ii,jj]*3 / twostar_mb[ii,jj]
      print(c(kk,ii,jj,"SS"))
    }
    ed=Sys.time()
    tiktok_sub[kk,ii]=ed-st
    
   
    tri_bt_sub[kk,ii,]=tri_mb[ii,]
    ts_bt_sub[kk,ii,]=twostar_mb[ii,]
    trans_bt_sub[kk,ii,]=trans_mb[ii,]
    
    ################### EG ##################################
    st=Sys.time()
    for(jj in 1:B){
      idx=sample(n,n,replace = T)
      P=A[idx,idx]
      P2=eigenMatMult(P,P)
      diag(P2)=0
      P3= P2* P
      num_tri_i2=colSums(P3)/2
      
      tri_mb[ii,jj]= (sum(num_tri_i2)/3)/(choose(n,3)*c_hat^3*tau_n^3)
      
      twostar_tot2=sum(colSums(P2*(1-P)))/2
      twostar_mb[ii,jj]= twostar_tot2/(choose(n,3)*c_hat^2*tau_n^2*(1-c_hat*tau_n))
      
      trans_mb[ii,jj]=tri_mb[ii,jj] *3 / twostar_mb[ii,jj]
      print(c(kk,ii,jj,'EG'))
    }
    ed=Sys.time()
    tiktok_eg[kk,ii]=ed-st
    
    
    tri_bt_eg[kk,ii,]=tri_mb[ii,]
    ts_bt_eg[kk,ii,]=twostar_mb[ii,]
    trans_bt_eg[kk,ii,]=trans_mb[ii,]
    
    
    #######################################################
    #fast-sketched
    
    ############# using our ustatistic sampling method #########
    N=floor(50*log(n))
    r=3
    st=Sys.time()
    M=cal_tri_h1(n,N,A) # see source code function
    h1=M[,N]/(c_hat^3*tau_n^3)
    
    
    Tnhat_us[kk,ii]=mean(h1)
    g1=h1-Tnhat_us[kk,ii]
    Sn_hat_formb_us= sqrt(r^2/n^2*sum(g1^2))
    
    
    Mts=cal_vstar_h1(n,N,A)
    h1_ts=Mts[,N]/(c_hat^2*tau_n^2*(1-c_hat*tau_n))
    
    Tnhat_us_ts[kk,ii]=mean(h1_ts)
    g1_ts=h1_ts-Tnhat_us_ts[kk,ii]
    Sn_hat_formb_us_ts= sqrt(r^2/n^2*sum(g1_ts^2))
    
    Tnhat_us_trans[kk,ii]=Tnhat_us[kk,ii]*3/Tnhat_us_ts[kk,ii]
    
    ed=Sys.time()
    time_us[kk,ii]=ed-st
    
    ############## using Kato method ####################
    N3=floor(n*50*log(n)/r)
    st=Sys.time()
    h1_kato=cal_kato_tri_h1(n,N3,A)/(N3*c_hat^3*tau_n^3)#-Tnhat_kato[ii,tt]
    # see source code function
    
    
    Tnhat_kato[kk,ii]=mean(h1_kato)
    g1kato=h1_kato-Tnhat_kato[kk,ii]
    Sn_hat_formb_kato= sqrt(r^2/n^2*sum(g1kato^2))
    
    h1_kato_ts=cal_kato_vstar_h1(n,N3,A)/(N3*c_hat^2*tau_n^2*(1-c_hat*tau_n))
    Tnhat_kato_ts[kk,ii]=mean(h1_kato_ts)
    g1kato_ts=h1_kato_ts-Tnhat_kato_ts[kk,ii]
    Sn_hat_formb_kato_ts= sqrt(r^2/n^2*sum(g1kato_ts^2))
    
    Tnhat_kato_trans[kk,ii]=Tnhat_kato[kk,ii]*3/Tnhat_kato_ts[kk,ii]
    
    ed=Sys.time()
    time_kato[kk,ii]=ed-st
    
    
    #MB-L-apx
    st=Sys.time()
    for(jj in 1:B){
      ee1= rnorm(n,1,sqrt(1/2)) 
      ee2= rnorm(n,1,sqrt(1/3))
      e_vec=ee1 * ee2
      TnLboot_us[ii,jj]=Tnhat_us[kk,ii]+ 3/n*sum((e_vec-1)*g1)
      TnLboot_kato[ii,jj]=Tnhat_kato[kk,ii]+ 3/n*sum((e_vec-1)*g1kato)
      
      TnLboot_us_ts[ii,jj]=Tnhat_us_ts[kk,ii]+ 3/n*sum((e_vec-1)*g1_ts)
      TnLboot_kato_ts[ii,jj]=Tnhat_kato_ts[kk,ii]+ 3/n*sum((e_vec-1)*g1kato_ts)
      
      TnLboot_us_trans[ii,jj]=3*TnLboot_us[ii,jj]/TnLboot_us_ts[ii,jj]
      TnLboot_kato_trans[ii,jj]=3*TnLboot_kato[ii,jj]/TnLboot_kato_ts[ii,jj]
      
      print(c(kk,ii,jj))
    }
    ed=Sys.time()
    tiktok_btl[kk,ii]=ed-st
    
    tri_bt_us[kk,ii,]= TnLboot_us[ii,]
    tri_bt_kato[kk,ii,]= TnLboot_kato[ii,]
    
    ts_bt_us[kk,ii,]=TnLboot_us_ts[ii,]
    ts_bt_kato[kk,ii,]=TnLboot_kato_ts[ii,]
    
    trans_bt_us[kk,ii,]=TnLboot_us_trans[ii,]
    trans_bt_kato[kk,ii,]=TnLboot_kato_trans[ii,]
     
    res=list(
      dg_hat=dg_hat,
      tri_hat=tri_hat,
      twostar_hat=twostar_hat,
      trans_hat=trans_hat,
    
      # time
      tiktok_pre=tiktok_pre,
      tiktok_ew=tiktok_ew,
      tiktok_gauprod=tiktok_gauprod,
      tiktok_linqua=tiktok_linqua,
      tiktok_add=tiktok_add,
      tiktok_lz=tiktok_lz,
      tiktok_sub=tiktok_sub,
      tiktok_eg=tiktok_eg,
      
      tiktok_btl=tiktok_btl,
      
      time_us=time_us,
      time_kato=time_kato,
      Tnhat_us=Tnhat_us,
      Tnhat_kato=Tnhat_kato,
      Tnhat_us_ts=Tnhat_us_ts,
      Tnhat_kato_ts=Tnhat_kato_ts,
      Tnhat_us_trans=Tnhat_us_trans,
      Tnhat_kato_trans=Tnhat_kato_trans,
      
      #### save bootstrap samples
      tri_bt_ew=tri_bt_ew,
      tri_bt_linqua= tri_bt_linqua,
      tri_bt_add=tri_bt_add,
      tri_bt_lz=tri_bt_lz,
      tri_bt_sub=tri_bt_sub,
      tri_bt_eg=tri_bt_eg,
      tri_bt_us=tri_bt_us,
      tri_bt_kato=tri_bt_kato,
      
      ts_bt_ew=ts_bt_ew,
      ts_bt_linqua= ts_bt_linqua,
      ts_bt_add=ts_bt_add,
      ts_bt_lz=ts_bt_lz,
      ts_bt_sub=ts_bt_sub,
      ts_bt_eg=ts_bt_eg,
      ts_bt_us=ts_bt_us,
      ts_bt_kato=ts_bt_kato,
      
      trans_bt_ew=trans_bt_ew,
      trans_bt_linqua= trans_bt_linqua,
      trans_bt_add=trans_bt_add,
      trans_bt_lz=trans_bt_lz,
      trans_bt_sub=trans_bt_sub,
      trans_bt_eg=trans_bt_eg,
      trans_bt_us=trans_bt_us,
      trans_bt_kato=trans_bt_kato,
      
      #### save corrections res
      tri_p1mat=tri_p1mat,
      tri_q1mat=tri_q1mat,
      tri_crtmat=tri_crtmat,
      tri_Eg13mat=tri_Eg13mat,
      tri_Eg1g1g2mat=tri_Eg1g1g2mat,
      tri_xi1mat=tri_xi1mat,
      
      ts_p1mat=ts_p1mat,
      ts_q1mat=ts_q1mat,
      ts_crtmat=ts_crtmat,
      ts_Eg13mat=ts_Eg13mat,
      ts_Eg1g1g2mat=ts_Eg1g1g2mat,
      ts_xi1mat=ts_xi1mat,
      
     
      p1trans_mat=p1trans_mat,
      q1trans_mat=q1trans_mat,
      trans_crtmat=trans_crtmat
      
    )
    
    if(input==1){
      save(res,file="coverage-sbm-6taus-all-n500.RData")
    }
    if(input==2){
      save(res,file="coverage-smg-6taus-all-n500.RData")
    }
  }
  
}




