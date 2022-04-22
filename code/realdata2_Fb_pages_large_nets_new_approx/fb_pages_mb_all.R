rm(list=ls())
set.seed(123)

setwd("../source/")

Rcpp::sourceCpp("find_neighbors.cpp")
Rcpp::sourceCpp("brute_counts_4cyc_3pth.cpp")
Rcpp::sourceCpp("approx_counts_us_newsp.cpp")

setwd("../realdata2_Fb_pages_large_nets_new_approx/")

library(igraph)
library(Matrix)
library(spam)

########### read in data ##############################################
dat_ath=read.csv("facebook_clean_data/athletes_edges.csv")

dat_pf=read.csv("facebook_clean_data/public_figure_edges.csv")

dat_pol=read.csv("facebook_clean_data/politician_edges.csv")

dat_art=read.csv("facebook_clean_data/artist_edges.csv")

mat_list_orig=list(dat_ath,dat_pf,dat_pol,dat_art)
n_vec=vector()

graph_list=list()
mat_list=list()
for(i in 1:4){
  dat1=as.matrix(mat_list_orig[[i]])
  colnames(dat1)=NULL
  graph_list[[i]]=graph_from_edgelist(dat1+1,directed = F)
  mat1=as_adjacency_matrix(graph_list[[i]],sparse = T)
  diag(mat1)=0
  mat1[mat1==2]=1
  mat_list[[i]]=mat1
  n_vec[i]=dim(mat1)[1]
}

rhon_hat=vector()
for(i in 1:4){
  n=n_vec[i]
  A=mat_list[[i]]
  rhon_hat[i]=(sum(colSums(A))/2)/(choose(n,2))
}

############### create space holders ##################################

B=1000
Tnhat_ut_tri=vector()
Tnhat_ut_ts=vector()
Tnhat_ut_fcq=vector()
Tnhat_ut_pth=vector()

timeg1s_tri=matrix(0,1,4)
timeg1s_ts=matrix(0,1,4)
timeg1s_fcq=matrix(0,1,4)
timeg1s_pth=matrix(0,1,4)

TnLboot_ut_tri=matrix(0,4,B)
TnLboot_ut_ts=matrix(0,4,B)
TnLboot_ut_fcq=matrix(0,4,B)
TnLboot_ut_pth=matrix(0,4,B)

tiktok_btl = vector()


mbvec_ut_tri=matrix(0,4,B)
mbvec_ut_ts=matrix(0,4,B)
mbvec_ut_fcq=matrix(0,4,B)
mbvec_ut_pth=matrix(0,4,B)

var_ut_tri=vector()
var_ut_ts=matrix(0,4,B)
var_ut_fcq=matrix(0,4,B)
var_ut_pth=matrix(0,4,B)

timeg1s_tri_old=matrix(0,1,4)
timeg1s_ts_old=matrix(0,1,4)
timeg1s_fcq_old=matrix(0,1,4)
timeg1s_pth_old=matrix(0,1,4)


#############  i is network index, there are four networks here ####
for(i in 1:4){
  # for each network i
  print(i)
  
  # from read in data
  n=n_vec[i]
  A=mat_list[[i]]
  n=dim(A)[1]
  dg=colSums(A)
  rhon=(sum(dg)/2)/(choose(n,2))
  
  # set N 
  
  N=floor(50*log(n))
  
  # find neighbours using cpp function 
  nbs_list=find_all(n,A)
  
  ###
  st=Sys.time()
  M=cal_3path_h1_newsp(n,N,A,nbs_list,dg)
  g1_pth_uncentered=M/(rhon^3)
  ed=Sys.time()
  timeg1s_pth[1,i]=difftime(ed, st, units = "secs")
  
  r=4
  Tnhat_ut_pth[i]=mean(g1_pth_uncentered)
  g1_pth=g1_pth_uncentered - Tnhat_ut_pth[i]
  
  ### four cliques
  # st=Sys.time()
  # M=cal_4clique_h1_newsp(n,N,A,nbs_list,dg)
  # g1_fcq_uncentered=M/(rhon^6)
  # ed=Sys.time()
  # timeg1s_fcq[1,i]=ed-st
  # 
  # r=4
  # Tnhat_ut_fcq[i]=mean(g1_fcq_uncentered)
  # g1_fcq=g1_fcq_uncentered-Tnhat_ut_fcq[i]
  # Sn_hat_formb_ut_fcq= sqrt(r^2/n^2*sum(g1_fcq^2))


  ### triangles
  st=Sys.time()
  M=cal_tri_h1_newsp(n,N,A,nbs_list,dg)
  g1_tri_uncentered=M/(rhon^3)
  ed=Sys.time()
  timeg1s_tri[1,i]=difftime(ed, st, units = "secs")

  r=3
  Tnhat_ut_tri[i]=mean(g1_tri_uncentered)
  g1_tri=g1_tri_uncentered-Tnhat_ut_tri[i]
  
  
  ### v-stars
  st=Sys.time()
  M=cal_vstar_h1_newsp(n,N,A,nbs_list,dg)
  g1_ts_uncentered=M/(rhon^2)
  ed=Sys.time()
  timeg1s_ts[1,i]=difftime(ed, st, units = "secs")
  
  r=3
  Tnhat_ut_ts[i]=mean(g1_ts_uncentered)
  g1_ts=g1_ts_uncentered-Tnhat_ut_ts[i]
  
  ######################################### if want to compare to old 
  ### old method time
  
  # # old 3 path
  # #old_count_3path=count_3path(n,as.matrix(A))/(choose(n,4)*rhon^3)
  # st=Sys.time()
  # old_count_3path_g1=count_3path_i(n,A)/(choose(n,4)*rhon^3)
  # old_count_3path=sum(old_count_3path_g1)/4
  # ed=Sys.time()
  # timeg1s_pth_old[1,i]=difftime(ed, st, units = "secs")
  # 
  # 
  # # old 4 cliques
  # #old_4c_count=count_4cliques(n,as.matrix(A))/(choose(n,4)*rhon^6)
  # 
  # #old triangles
  # st=Sys.time()
  # A2= A %*% A
  # diag(A2)=0
  # A3= A2* 1*(A==1)
  # num_tri_temp=colSums(A3)/2
  # 
  # (sum(num_tri_temp)/3)/(choose(n,3)*rhon^3)
  # ed=Sys.time()
  # timeg1s_tri_old[1,i]=difftime(ed, st, units = "secs")
  # 
  # # old v-star
  # st=Sys.time()
  # A2= A %*% A
  # diag(A2)=0
  # (sum(colSums(A2*(1-A)))/2)/(choose(n,3)*rhon^2)
  # ed=Sys.time()
  # timeg1s_ts_old[1,i]=difftime(ed, st, units = "secs")
  # #############################################
  # 
  
  
  #MB-L-apx
  st=Sys.time()
  for(jj in 1:B){
    ee1= rnorm(n,1,sqrt(1/2))
    ee2= rnorm(n,1,sqrt(1/3))
    e_vec=ee1 * ee2
    ### note here N is not >> nrho^s, so here we do 1 instead of r in multiplier
    TnLboot_ut_pth[i,jj]=Tnhat_ut_pth[i]+ 1/n*sum((e_vec-1)*g1_pth)
    #TnLboot_ut_fcq[i,jj]=Tnhat_ut_fcq[i]+ 4/n*sum((e_vec-1)*g1_fcq)
    TnLboot_ut_tri[i,jj]=Tnhat_ut_tri[i]+ 1/n*sum((e_vec-1)*g1_tri)
    TnLboot_ut_ts[i,jj]=Tnhat_ut_ts[i]+ 1/n*sum((e_vec-1)*g1_ts)
    print(c(i,jj))
  }
  ed=Sys.time()
  tiktok_btl[i]=ed-st
  
  mbvec_ut_pth[i,]=  TnLboot_ut_pth[i,]
  #mbvec_ut_fcq[i,]=  TnLboot_ut_fcq[i,]
  mbvec_ut_tri[i,]=  TnLboot_ut_tri[i,]
  mbvec_ut_ts[i,]=  TnLboot_ut_ts[i,]
  
  var_ut_pth[i]=var(TnLboot_ut_pth[i,])
  #var_ut_fcq[i]=var(TnLboot_ut_fcq[i,])
  var_ut_fcq[i]=var(TnLboot_ut_tri[i,])
  var_ut_ts[i]=var(TnLboot_ut_ts[i,])
  
}



res=list(timeg1s_pth=timeg1s_pth,
         timeg1s_tri=timeg1s_tri,
         timeg1s_ts=timeg1s_ts,
         mbvec_ut_pth=mbvec_ut_pth,
         mbvec_ut_tri=mbvec_ut_tri,
         mbvec_ut_ts=mbvec_ut_ts,
         var_ut_pth=var_ut_pth,
         var_ut_tri=var_ut_tri,
         var_ut_ts=var_ut_ts,
         Tnhat_ut_pth=Tnhat_ut_pth,
         Tnhat_ut_tri=Tnhat_ut_tri,
         Tnhat_ut_ts=Tnhat_ut_ts,
         n_vec=n_vec,
         rhon_hat=rhon_hat)


rm(mat_list)


save(res,file='fb-pages-all-mblapx-new.RData')


timecal <- function(n,rho){
  return(n*50*log(n)*n*rho*log(n*rho))
}


btci_ut_pth=matrix(0,4,2)
for(i in 1:4){
  btci_ut_pth[i,]=c(quantile(res$mbvec_ut_pth[i,],0.05/3),quantile(res$mbvec_ut_pth[i,],1-0.05/3))
}

pdf("BonfCo_FB_CI_3grp_3path_new.pdf",width=5, height = 3)
par(mar=c(3,3,1,1))
plot(c(1,2,3),res$Tnhat_ut_pth[c(1,3,4)],xlim=c(0,4),ylim=c(10,250),pch=20,xlab="",ylab="",cex=0.5,xaxt="n",cex.axis=1)
points(c(1,2,3),btci_ut_pth[c(1,3,4),1],pch="-",cex=1)
points(c(1,2,3),btci_ut_pth[c(1,3,4),2],pch="-",cex=1)
segments(1,btci_ut_pth[1,1],1,btci_ut_pth[1,2],lty=1)
#segments(2,btci_ut_pth[2,1],2,btci_ut_pth[2,2],lty=1)
segments(2,btci_ut_pth[3,1],2,btci_ut_pth[3,2],lty=1)
segments(3,btci_ut_pth[4,1],3,btci_ut_pth[4,2],lty=1)

axis(1,cex.axis=1, at=1, labels = 'Athletes')
#axis(1,cex.axis=0.6, at=2, labels = 'Public Figures')
axis(1,cex.axis=1, at=2, labels = 'Politicians')
axis(1,cex.axis=1, at=3, labels = 'Artists')

dev.off()

btci_ut_tri=matrix(0,4,2)
for(i in 1:4){
  btci_ut_tri[i,]=c(quantile(res$mbvec_ut_tri[i,],0.05/3),quantile(res$mbvec_ut_tri[i,],1-0.05/3))
}

pdf("BonfCo_FB_CI_3grp_tri_new.pdf",width=5, height = 3)
par(mar=c(2,2,1,1))
plot(c(1,2,3),res$Tnhat_ut_tri[c(1,3,4)],xlim=c(0,4),ylim=c(200,600),pch=20,xlab="",ylab="",cex=0.5,xaxt="n",cex.axis=1)
points(c(1,2,3),btci_ut_tri[c(1,3,4),1],pch="-",cex=1)
points(c(1,2,3),btci_ut_tri[c(1,3,4),2],pch="-",cex=1)
segments(1,btci_ut_tri[1,1],1,btci_ut_tri[1,2],lty=1)
#segments(2,btci_ut_tri[2,1],2,btci_ut_tri[2,2],lty=1)
segments(2,btci_ut_tri[3,1],2,btci_ut_tri[3,2],lty=1)
segments(3,btci_ut_tri[4,1],3,btci_ut_tri[4,2],lty=1)

axis(1,cex.axis=1, at=1, labels = 'Athletes')
#axis(1,cex.axis=0.6, at=2, labels = 'Public Figures')
axis(1,cex.axis=1, at=2, labels = 'Politicians')
axis(1,cex.axis=1, at=3, labels = 'Artists')

dev.off()

btci_ut_ts=matrix(0,4,2)
for(i in 1:4){
  btci_ut_ts[i,]=c(quantile(res$mbvec_ut_ts[i,],0.05/3),quantile(res$mbvec_ut_ts[i,],1-0.05/3))
}

pdf("BonfCo_FB_CI_3grp_vstar_new.pdf",width=5, height = 3)
par(mar=c(3,3,1,1))
plot(c(1,2,3),res$Tnhat_ut_ts[c(1,3,4)],xlim=c(0,4),ylim=c(0,20),pch=20,xlab="",ylab="",cex=0.5,xaxt="n",cex.axis=1)
points(c(1,2,3),btci_ut_ts[c(1,3,4),1],pch="-",cex=1)
points(c(1,2,3),btci_ut_ts[c(1,3,4),2],pch="-",cex=1)
segments(1,btci_ut_ts[1,1],1,btci_ut_ts[1,2],lty=1)
#segments(2,btci_ut_ts[2,1],2,btci_ut_ts[2,2],lty=1)
segments(2,btci_ut_ts[3,1],2,btci_ut_ts[3,2],lty=1)
segments(3,btci_ut_ts[4,1],3,btci_ut_ts[4,2],lty=1)

axis(1,cex.axis=1, at=1, labels = 'Athletes')
#axis(1,cex.axis=0.6, at=2, labels = 'Public Figures')
axis(1,cex.axis=1, at=2, labels = 'Politicians')
axis(1,cex.axis=1, at=3, labels = 'Artists')
dev.off()

