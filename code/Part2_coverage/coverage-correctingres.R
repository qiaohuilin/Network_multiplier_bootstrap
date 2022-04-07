setwd("../Part2_coverage/")

##############load result from the coverage-truth.R and coverage-sim.R#####
#load truth
load("coverage-truth-6taus.RData")
#load result from coverage-sim.R
load("coverage-sbm-6taus-all-n500.RData")

covg_sbm_triangles=matrix(0,6,6)

# corrected 
for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,1]=covg_sbm_triangles[kk,1] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_linqua[kk,ii,],0.025)+res$tri_crtmat[kk,ii] 
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_linqua[kk,ii,],0.975)+res$tri_crtmat[kk,ii])
    
  }
}


for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,4]=covg_sbm_triangles[kk,4] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_sub[kk,ii,],0.025)+res$tri_crtmat[kk,ii] 
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_sub[kk,ii,],0.975)+res$tri_crtmat[kk,ii])
    
  }
}

for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,5]=covg_sbm_triangles[kk,5] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_eg[kk,ii,],0.025)+res$tri_crtmat[kk,ii] 
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_eg[kk,ii,],0.975)+res$tri_crtmat[kk,ii])
    
  }
}

# methods that are not higher order corrected
for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,2]=covg_sbm_triangles[kk,2] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_add[kk,ii,],0.025)
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_add[kk,ii,],0.975))
    
  }
}

for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,3]=covg_sbm_triangles[kk,3] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_lz[kk,ii,],0.025)
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_lz[kk,ii,],0.975))
    
  }
}

for(kk in 1:6){
  for(ii in 1:200){
    covg_sbm_triangles[kk,6]=covg_sbm_triangles[kk,6] + 
      1*(true$numtri_hat_true[1,kk] > quantile(res$tri_bt_ut[kk,ii,],0.025) 
         & true$numtri_hat_true[1,kk] < quantile(res$tri_bt_ut[kk,ii,],0.975))
    
  }
}

# repeat same for v-star and transitivity

covg_sbm_triangles=covg_sbm_triangles/2
covg_sbm_vs=covg_sbm_vs/2
covg_sbm_trans=covg_sbm_trans/2

# same for smooth graphon


res=list(covg_sbm_triangles=covg_sbm_triangles,covg_sbm_triangles_old=covg_sbm_triangles_old,covg_sbm_vs=covg_sbm_vs,covg_sbm_vs_old=covg_sbm_vs_old,
         covg_sbm_trans=covg_sbm_trans,covg_sbm_trans_old=covg_sbm_trans_old
)



