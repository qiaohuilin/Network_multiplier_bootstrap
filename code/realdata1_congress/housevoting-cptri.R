############### note that in this example we only used MB-Q ########
############### these are small and dense graphs ###################
###################other methods are commented out ##################

setwd("../realdata1_congress/")

####################### run housevoting.R first ##################
#### load the A_list,rhon_hat,n_list, cped,cptri (all networks calculated there)
A_list=res$A_list
cped=res$cped
cptri=res$cptri

#spaceholders 
mbvec_linqua=array(0,c(32,B))
tiktok_pre=vector()
tiktok_linqua=vector()
var_linqua=vector()

B=2000
cped_mb=matrix(0,32,B)

#grids for cdf
x=seq(-3,3,by = 0.1)
spt=1:61

########### if you want to run other methods as well ##############
#load rcpp functions for approximate cross party counts (different from source)
#Rcpp::sourceCpp("../source/cross_party_approx_original.cpp")

# mbvec_ew=array(0,c(32,61))
# mbvec_gauprod=array(0,c(32,B))
# mbvec_linqua=array(0,c(32,B))
# mbvec_add=array(0,c(32,B))
# mbvec_lz=array(0,c(32,B))
# mbvec_sub=array(0,c(32,B))
# mbvec_eg=array(0,c(32,B))
# mbvec_us=array(0,c(32,B))
# mbvec_kato=array(0,c(32,B))
# 
# tiktok_pre=vector()
# tiktok_ew=vector()
# tiktok_gauprod=vector()
# tiktok_linqua=vector()
# tiktok_add=vector()
# tiktok_lz=vector()
# tiktok_sub=vector()
# tiktok_eg=vector()
# 
# 
# var_ew=vector()
# var_gauprod=vector()
# var_linqua=vector()
# var_add=vector()
# var_lz=vector()
# var_sub=vector()
# var_eg=vector()
# var_us=vector()
# var_kato=vector()
# Tnhat_us=vector()
# Tnhat_kato=vector()
# timeg1=matrix(0,2,32) #first row us, second row kato
# tiktok_btl = vector()
# 
# TnLboot_us=matrix(0,32,B)
# TnLboot_kato=matrix(0,32,B)
# 

################# if use linear methods the corrections are here #######
# cped_p1mat=vector()
# cped_q1mat=vector()
# cped_crtmat=vector()

# #p1 and q1 are symmetric re x, z_alpha and z_1-alpha have same p1 q1 value
# p1 <- function(x,xi1,r,Eg13,Eg1g1g2){
#   return(-(x^2-1)/(6*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2))
# }
# 
# 
# q1 <- function(x,xi1,r,Eg13,Eg1g1g2){
#   return(1/(xi1^3)*((2*x^2+1)/6*Eg13+(r-1)*(x^2+1)/2*Eg1g1g2))
# }

################functions for edgeworth###############################

# EWCDF_sd <- function(x,n,xi1,r,Eg13,Eg1g1g2){
#   gx=pnorm(x)-dnorm(x)*(x^2-1)/(6*sqrt(n)*(xi1^3))*(Eg13+3*(r-1)*Eg1g1g2)
#   return(gx)
# }

##############################################################3

for(i in 1:32){
  A=A_list[[i]]
  n=dim(A)[1]
  rhon=rhon_hat[i]
  ###### you can load cptri from previous result 
  ###### or if you want to recalculate it here 
  
  # A2=eigenMatMult(A,A)
  # A3= A2* 1*(A==1)
  # num_tri_temp=colSums(A3)/2
  # num_tri_tot=sum(num_tri_temp)/3
  # Asp1=A[party_label[[i]]==100,party_label[[i]]==100]
  # Asp2=A[party_label[[i]]==200,party_label[[i]]==200]
  # A3sp1=eigenMatMult(Asp1,Asp1)*Asp1
  # A3sp2=eigenMatMult(Asp2,Asp2)*Asp2
  # numtri_sp1=colSums(A3sp1)/2
  # numtri_sp2=colSums(A3sp2)/2
  # numtri_cp <- vector()
  # numtri_cp[party_label[[i]]==100]=num_tri_temp[party_label[[i]]==100]-numtri_sp1
  # numtri_cp[party_label[[i]]==200]=num_tri_temp[party_label[[i]]==200]-numtri_sp2
  # numtri_cp_tot=sum(numtri_cp)/3
  # cptri[i]=numtri_cp_tot/(choose(n,3)*rhon^3)
  
  ################## pre-calculations #######################
  st=Sys.time()
  r=3
  g1=numtri_cp/(choose(n-1, 2)*rhon^3)-cptri[i]
  Sn_hat_formb= sqrt(r^2/n^2*sum(g1^2))
  h2mat=A3
  h2mat[party_label[[i]]==100,party_label[[i]]==100]=g2mat[party_label[[i]]==100,party_label[[i]]==100]-A3sp1
  h2mat[party_label[[i]]==200,party_label[[i]]==200]=g2mat[party_label[[i]]==200,party_label[[i]]==200]-A3sp2
  g2=h2mat/(choose(n-2,1)*rhon^3)-cptri[i]-g1*matrix(1,n,n)- t(g1*matrix(1,n,n))
  diag(g2)=0
  g2tilde=g2mat/(choose(n-2,1)*rhon^3)-cptri[i]
  diag(g2tilde)=0
  
  ed=Sys.time()
  tiktok_pre[i]=ed-st
  
  ###################### MB-Q #################################
  st=Sys.time()
  for(jj in 1:B){
    ee1= rnorm(n,1,sqrt(1/2)) 
    ee2= rnorm(n,1,sqrt(1/3))
    e_vec=ee1 * ee2
    
    wg2mat=(e_vec*matrix(1,n,n) * t(e_vec*matrix(1,n,n)) - matrix(1,n,n)) * g2tilde - (e_vec*g1*matrix(1,n,n)) - t(e_vec*g1*matrix(1,n,n)) 
    diag(wg2mat)=0
    wg2=sum(wg2mat)/2
    cptri_mb[i,jj]=cptri[i]+r/n*sum((e_vec-1)*g1)+r*(r-1)/(n*(n-1))*wg2 
    print(c(i,jj))
  }
  ed=Sys.time()
  tiktok_linqua[i]=ed-st
  
  Tmb_lq=(cptri_mb[i,]-mean(cptri_mb[i,]))/Sn_hat_formb #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
  Tmbcdf_lq=sapply(1:length(x),function(qq) mean(Tmb_lq<x[qq]))
  mbvec_linqua[i,]= cptri_mb[i,]
  var_linqua[i]=var(cptri_mb[i,])
  
  ###################if you want to use other methods ##########################
  
  #################### Edgeworth Expansion ##########################
  # st=Sys.time()
  # Eg13= mean(g1^3)
  # Eg1g1g2=sum(g2*(g1*matrix(1,n,n))* t(g1*matrix(1,n,n)))/(2*choose(n,2))
  # xi1=sqrt(mean(g1^2))
  # 
  # ew=EWCDF_sd(x,n,xi1,3,Eg13,Eg1g1g2)
  # cptri_p1mat[i]=p1(1.96,xi1,3,Eg13,Eg1g1g2)
  # cptri_q1mat[i]=q1(1.96,xi1,3,Eg13,Eg1g1g2)
  # cptri_crtmat[i]=n^(-1)*Sn_hat_formb*( cptri_p1mat[i]+ cptri_q1mat[i])
  # 
  # 
  # mbvec_ew[i,]=ew[spt]
  # 
  # ed=Sys.time()
  # tiktok_ew[i]=ed-st
  
  
  ###################### MB-L ################################
  # st=Sys.time()
  # for(jj in 1:B){
  #   ee1= rnorm(n,1,sqrt(1/2)) 
  #   ee2= rnorm(n,1,sqrt(1/3))
  #   e_vec=ee1 * ee2
  #   cptri_mb[i,jj]=cptri[i]+r/n*sum((e_vec-1)*g1)
  #   print(c(i,jj))
  # }
  # ed=Sys.time()
  # tiktok_add[i]=ed-st
  # 
  # Tmb_add=(cptri_mb[i,]-mean(cptri_mb[i,]))/Sn_hat_formb #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
  # Tmbcdf_add=sapply(1:length(x),function(qq) mean(Tmb_add<x[qq]))
  # mbvec_add[i,]= cptri_mb[i,]
  # var_add[i]=var(cptri_mb[i,])
  
  ###################### LS ##################################
  # st=Sys.time()
  # eg=eigen(A)
  # idx=order(abs(eg$values),decreasing = TRUE) 
  # qct=abs(eg$values) > 2.01*sqrt(n*mean(A))
  # egval=eg$values[qct]
  # egvec=eg$vectors[,qct]
  # if(sum(qct)==1){
  #   egvec=matrix(egvec, n, 1)
  #   P= (egvec * egval)  %*% t(egvec)
  # }else{
  #   P= egvec %*% diag(egval) %*% t(egvec)
  # }
  # 
  # diag(P)=0
  # P[P>1]=1
  # P[P<0]=0
  # 
  # P2=eigenMatMult(P,P)
  # diag(P2)=0
  # P3= P2*P
  # num_tri_temp=colSums(P3)/2
  # num_tri_tot=sum(num_tri_temp)/3
  # Psp1=P[party_label[[i]]==100,party_label[[i]]==100]
  # Psp2=P[party_label[[i]]==200,party_label[[i]]==200]
  # P3sp1=eigenMatMult(Psp1,Psp1)*Psp1
  # P3sp2=eigenMatMult(Psp2,Psp2)*Psp2
  # numtri_sp1=colSums(P3sp1)/2
  # numtri_sp2=colSums(P3sp2)/2
  # 
  # numtri_cp <- vector()
  # numtri_cp[party_label[[i]]==100]=
  #   num_tri_temp[party_label[[i]]==100]-numtri_sp1
  # numtri_cp[party_label[[i]]==200]=
  #   num_tri_temp[party_label[[i]]==200]-numtri_sp2
  # 
  # numtri_cp_tot=sum(numtri_cp)/3
  # 
  # cptrilz= numtri_cp_tot/(choose(n,3)*rhon^3)
  
  # cptrilz_i=numtri_cp/(choose(n-1,2)*rhon^3)
  # g1lz=cptrilz_i-cptrilz
  # Sn_hat_lz= sqrt(r^2/n^2*sum(g1lz^2))
  # 
  # W= rmultinom(B,size=n,prob = rep(1/n,n))
  # for(jj in 1: B) {
  #   cptri_mb[i,jj]=  cptrilz+ 3*(sum(W[,jj]*cptrilz_i)/n- cptrilz) #3*(sum(W[,jj]*tri_i)/n) # #(sum(diag(A3))/6)/(choose(n,3)*rhon^3)
  #   print(c(i,jj))
  # }
  # ed=Sys.time()
  # tiktok_lz[i]=ed-st
  # 
  # Tmb_lz=(cptri_mb[i,]-mean(cptri_mb[i,]))/Sn_hat_lz #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
  # Tmbcdf_lz=sapply(1:length(x),function(qq) mean(Tmb_lz<x[qq]))
  # mbvec_lz[i,]= cptri_mb[i,]
  # var_lz[i]=var(cptri_mb[i,])
  
  ####################### SS #############################
  # st=Sys.time()
  # Sn_hat_formbsb=NULL
  # for(jj in 1:B){
  #   b=floor(0.5*n)
  #   idx=sample(1:n,b,replace = F)
  #   subA=A[idx,idx]
  #   sublabel=party_label[[i]][idx]
  #   
  #   subA2=eigenMatMult(subA,subA)
  #   subA3= subA2* 1*(subA==1)
  #   num_tri_temp=colSums(subA3)/2
  #   num_tri_tot=sum(num_tri_temp)/3
  #   subAsp1=subA[sublabel==100,sublabel==100]
  #   subAsp2=subA[sublabel==200,sublabel==200]
  #   subA3sp1=eigenMatMult(subAsp1,subAsp1)*subAsp1
  #   subA3sp2=eigenMatMult(subAsp2,subAsp2)*subAsp2
  #   numtri_sp1=colSums(subA3sp1)/2
  #   numtri_sp2=colSums(subA3sp2)/2
  #   numtri_cp <- vector()
  #   numtri_cp[sublabel==100]=num_tri_temp[sublabel==100]-numtri_sp1
  #   numtri_cp[sublabel==200]=num_tri_temp[sublabel==200]-numtri_sp2
  #   numtri_cp_tot=sum(numtri_cp)/3
  #   
  #   
  #   cptri_mb[i,jj]=numtri_cp_tot/(choose(b,3)*rhon^3)
  #   
  #   g1sb=numtri_cp/(choose(b-1,2)*rhon^3)-cptri_mb[i,jj]
  #   Sn_hat_formbsb= sqrt(r^2/b^2*sum(g1sb^2))
  #   print(c(i,jj))
  # }
  # ed=Sys.time()
  # tiktok_sub[i]=ed-st
  # 
  # Tmb_sub=(cptri_mb[i,]-mean(cptri_mb[i,]))/mean(Sn_hat_formbsb*sqrt(b/n))  #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
  # Tmbcdf_sub=sapply(1:length(x),function(qq) mean(Tmb_sub<x[qq]))
  # mbvec_sub[i,]= cptri_mb[i,]
  # var_sub[i]=var(cptri_mb[i,])*b/n
  
  ########################### EG ################################
  # st=Sys.time()
  # for(jj in 1:B){
  #   idx=sample(n,n,replace = T)
  #   P=A[idx,idx]
  #   labeleg=party_label[[i]][idx]
  #   
  #   P2=eigenMatMult(P,P)
  #   diag(P2)=0
  #   P3=P2*P
  #   
  #   num_tri_temp=colSums(P3)/2
  #   num_tri_tot=sum(num_tri_temp)/3
  #   Psp1=P[labeleg==100,labeleg==100]
  #   Psp2=P[labeleg==200,labeleg==200]
  #   P3sp1=eigenMatMult(Psp1,Psp1)*Psp1
  #   P3sp2=eigenMatMult(Psp2,Psp2)*Psp2
  #   numtri_sp1=colSums(P3sp1)/2
  #   numtri_sp2=colSums(P3sp2)/2
  #   numtri_cp <- vector()
  #   numtri_cp[labeleg==100]=num_tri_temp[labeleg==100]-numtri_sp1
  #   numtri_cp[labeleg==200]=num_tri_temp[labeleg==200]-numtri_sp2
  #   numtri_cp_tot=sum(numtri_cp)/3
  #   
  #   cptri_mb[i,jj]=  numtri_cp_tot/(choose(n,3)*rhon^3)
  #   
  #   print(c(i,jj))
  # }
  # ed=Sys.time()
  # tiktok_eg[i]=ed-st
  # 
  # Tmb_eg=(cptri_mb[i,]-mean(cptri_mb[i,]))/Sn_hat_formb #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
  # Tmbcdf_eg=sapply(1:length(x),function(qq) mean(Tmb_eg<x[qq]))
  # mbvec_eg[i,]= cptri_mb[i,]
  # var_eg[i]=var(cptri_mb[i,])
  
  
  #######################################################
  #fast-sketched
  
#   N=floor(50*log(n))
#   r=3
#   st=Sys.time()
#   M=calcptrih1(n,N,A,party_label[[i]])
#   h1=M[,N]/(rhon^3)#-Tnhat[i,tt]
#   ed=Sys.time()
#   timeg1[1,i]=ed-st
#   
#   Tnhat_ut[i]=mean(h1)
#   g1=h1-Tnhat_ut[i]
#   Sn_hat_formb_ut= sqrt(r^2/n^2*sum(g1^2))
#   
#   k=3
#   N3=floor(n*50*log(n)/k)
#   
#   st=Sys.time()
#   h1kato=calcptrih1kato(n,N3,A,party_label[[i]])/(choose(n-1,2)*rhon^3)#-Tnhat_kato[i,tt]
#   ed=Sys.time()
#   timeg1[2,i]=ed-st
#   
#   Tnhat_kato[i]=mean(h1kato)
#   g1kato=h1kato-Tnhat_kato[i]
#   Sn_hat_formb_kato= sqrt(r^2/n^2*sum(g1kato^2))
#   
#   
#   
#   #MB-L-apx
#   st=Sys.time()
#   for(jj in 1:B){
#     ee1= rnorm(n,1,sqrt(1/2)) 
#     ee2= rnorm(n,1,sqrt(1/3))
#     e_vec=ee1 * ee2
#     TnLboot_ut[i,jj]=Tnhat_ut[i]+ 3/n*sum((e_vec-1)*g1)
#     TnLboot_kato[i,jj]=Tnhat_kato[i]+ 3/n*sum((e_vec-1)*g1kato)
#     print(c(i,jj))
#   }
#   ed=Sys.time()
#   tiktok_btl[i]=ed-st
#   
#   Tmb_ut=(TnLboot_ut[i,]-mean(TnLboot_ut[i,]))/Sn_hat_formb_ut #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
#   Tmbcdf_ut=sapply(1:length(x),function(qq) mean(Tmb_ut<x[qq]))
#   mbvec_ut[i,]=  TnLboot_ut[i,]
#   var_ut[i]=var(TnLboot_ut[i,])
#   
#   Tmb_kato=(TnLboot_kato[i,]-mean(TnLboot_kato[i,]))/Sn_hat_formb_kato #sqrt(var(cptri_mb[i,]))#Sn_hat_formb
#   Tmbcdf_kato=sapply(1:length(x),function(qq) mean(Tmb_kato<x[qq]))
#   mbvec_kato[i,]=  TnLboot_kato[i,]
#   var_kato[i]=var(TnLboot_kato[i,])

}

res=list(mbvec_linqua=mbvec_linqua,
         # mbvec_add=mbvec_add,
         # mbvec_ew=mbvec_ew,
         # mbvec_sub=mbvec_sub,
         # mbvec_lz=mbvec_lz,
         # mbvec_eg=mbvec_eg,
         # mbvec_us=mbvec_us,
         # mbvec_kato=mbvec_kato,
         
         var_linqua=var_linqua,
         # var_add=var_add,
         # var_ew=var_ew,
         # var_sub=var_sub,
         # var_lz=var_lz,
         # var_eg=var_eg,
         # var_us=var_us,
         # var_kato=var_kato,
         
         tiktok_pre=tiktok_pre,
         tiktok_linqua=tiktok_linqua,
         # tiktok_ew=tiktok_ew,
         # tiktok_add=tiktok_add,
         # tiktok_lz=tiktok_lz,
         # tiktok_sub=tiktok_sub,
         # tiktok_eg=tiktok_eg,
         # timeg1=timeg1,
         # tiktok_btl = tiktok_btl,
         # 
         # cped_p1mat=cped_p1mat,
         # cped_q1mat=cped_q1mat,
         # cped_crtmat=cped_crtmat,
         # 
         cptri=cptri
        )

save(res,file="MB-cptri-withcorrection.RData")

btci_linqua=matrix(0,32,2)
for(i in 1:32){
  #c(quantile(mbvec_add[i,],0.025)+cped_crtmat[i],quantile(mbvec_add[i,],0.975)+cped_crtmat[i])
  btci_linqua[i,]=c(quantile(mbvec_linqua[i,],0.05/32)+cptri_crtmat[i],quantile(mbvec_linqua[i,],1-0.05/32)+cptri_crtmat[i])
}
btci_add=matrix(0,32,2)
for(i in 1:32){
  #c(quantile(mbvec_add[i,],0.025),quantile(mbvec_add[i,],0.975))
  btci_add[i,]=c(quantile(mbvec_add[i,],0.05/32),quantile(mbvec_add[i,],1-0.05/32))
}



pdf("house-mb-cptri-bc.pdf",width=5,height=2.5)
par(mar=c(2,2,1,1))
plot(congress_list, res$cptri,type='l',ylab='CP Triangle Density',xlab='Congress',cex=0.5,ylim=c(0,2.3))
points(congress_list,res$cptri,pch=20,cex=0.7)
points(congress_list,res$cptri-3.16*sqrt(res$var_add),pch='-',cex=0.7)
points(congress_list,res$cptri+3.16*sqrt(res$var_add),pch='-',cex=0.7)
for(i in 1:32){
  segments(congress_list[i],res$cptri[i]-3.16*sqrt(res$var_add[i]),congress_list[i],res$cptri[i]+3.16*sqrt(res$var_add[i]))
}
dev.off()


pdf("house-mb-cptri-bc-quantile-mbadd.pdf",width=5,height=2.5)
par(mar=c(2,2,1,1))
plot(congress_list, rowMeans(btci_add),type='l',ylab='CP Edge Density',xlab='Congress',cex=0.5,ylim=c(0,2))
points(congress_list,rowMeans(btci_add),pch=20,cex=0.7)
points(congress_list,btci_add[,1],pch='-',cex=0.7)
points(congress_list,btci_add[,2],pch='-',cex=0.7)
for(i in 1:32){
  segments(congress_list[i],btci_add[i,1],congress_list[i],btci_add[i,2])
}
dev.off()

pdf("house-mb-cptri-bc-quantile-mblinqua.pdf",width=5,height=2.5)
par(mar=c(2,2,1,1))
plot(congress_list, rowMeans(btci_linqua),type='l',ylab='CP Edge Density',xlab='Congress',cex=0.5,ylim=c(0,2))
points(congress_list,rowMeans(btci_linqua),pch=20,cex=0.7)
points(congress_list,btci_linqua[,1],pch='-',cex=0.7)
points(congress_list,btci_linqua[,2],pch='-',cex=0.7)
for(i in 1:32){
  segments(congress_list[i],btci_linqua[i,1],congress_list[i],btci_linqua[i,2])
}
dev.off()

year_list=seq(1949,2011,by=2)
pdf("house-year-mb-cptri-bc-quantile-mblinqua.pdf",width=5,height=2.5)
par(mar=c(2,2,1,1))
plot(year_list, rowMeans(btci_linqua),type='l',ylab='CP Edge Density',xlab='Congress',cex=0.5,ylim=c(0,2),cex.axis=0.7)
points(year_list,rowMeans(btci_linqua),pch=20,cex=0.7)
points(year_list,btci_linqua[,1],pch='-',cex=0.7)
points(year_list,btci_linqua[,2],pch='-',cex=0.7)
for(i in 1:32){
  segments(year_list[i],btci_linqua[i,1],year_list[i],btci_linqua[i,2])
}
dev.off()
