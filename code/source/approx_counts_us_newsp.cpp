// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <vector>
#include <unordered_set>

using namespace Rcpp ;
using namespace arma ;




// [[Rcpp::export]]
arma::uvec get_sample_id_norep(int k, int n){
  arma::uvec id(k);
  arma::uvec ar(n);
  ar.fill(0);
  int ii;
  for (ii = 0; ii < k; ii++) {
    id(ii) = n * unif_rand();
    while(ar(id(ii))==1){
      id(ii) = n * unif_rand();
      //Rcout << "ii is: " <<  ii  << "; id(ii) is: " << id(ii) << " ; ar(id(ii)) is " << ar(id(ii)) << std::endl;;
      //cout id(ii) >> "ar(id(ii)):" >> ar(id(ii));
    }
    ar(id(ii))=1;
  }
  return id;
}



// [[Rcpp::export]]
arma::uvec get_sample_id_nov(int k, int n, IntegerVector v){
  arma::uvec id(k);
  int ii;
  for (ii = 0; ii < k; ii++) {
    id(ii) = n * unif_rand();
    while(v(id(ii))==1){
      id(ii) = n * unif_rand();
      //Rcout << "ii is: " <<  ii  << "; id(ii) is: " << id(ii) << " ; ar(id(ii)) is " << ar(id(ii)) << std::endl;;
      //cout id(ii) >> "ar(id(ii)):" >> ar(id(ii));
    }
  }
  return id;
}




// [[Rcpp::export]]
NumericVector cal_4clique_h1_newsp(int n, int n_sample, sp_mat A, List nbs_list, IntegerVector dg){
  NumericVector M(n);
  int k=4;
  int src=0;
  while(src<n){
    if(dg(src)<3){
      M(src)=0;
      src++;
    }else{
    double s=0;
    IntegerVector nbs = nbs_list[src];
    int nid=nbs.size();
    for(int j=1; j<n_sample; j++){
      arma::uvec S=get_sample_id_norep(nid,n-1);
      arma::uvec indexes=sort(S);
      arma::uvec sorted=sort_index(S);
      int i=0;
      while(i<(nid-2)){
        double t=0.0;
        int curr = indexes(i);
        int next = indexes(i+1);
        int next2 = indexes(i+2);
        int currnb= nbs(sorted(i));
        int nextnb=nbs(sorted(i+1));
        int nextnb2=nbs(sorted(i+2));
        if (curr % 3 ==1 && next==curr+1 && next2==curr+2){
          t=A(currnb,nextnb)*A(currnb,nextnb2)*A(nextnb,nextnb2);
          i=i+2;
        }
        s=s+t;
        i=i+1;
      }
    }
    s=s/((n-1)/(k-1));
    M(src)=s/n_sample;
    src++;
    }
    }
  return(M);
}




// [[Rcpp::export]]
NumericVector cal_tri_h1_newsp(int n, int n_sample, sp_mat A, List nbs_list, IntegerVector dg){
  NumericVector M(n);
  int k=3;
  int src=0;
  while(src<n){
    if(dg(src)<2){
      M(src)=0;
      src++;
    }else{
      double s=0;
      IntegerVector nbs = nbs_list[src];
      int nid=nbs.size();
      for(int j=1; j<n_sample; j++){
        arma::uvec S=get_sample_id_norep(nid,n-1);
        arma::uvec indexes=sort(S);
        arma::uvec sorted=sort_index(S);
        int i=0;
        while(i<(nid-1)){
          double t=0.0;
          int curr = indexes(i);
          int next = indexes(i+1);
          int currnb= nbs(sorted(i));
          int nextnb=nbs(sorted(i+1));
          if (curr % 2 ==1 && next==curr+1){
            t=A(currnb,nextnb);
            i=i+1;
          }
          s=s+t;
          i=i+1;
        }
      }
      s=s/((n-1)/(k-1));
      M(src)=s/n_sample;
      src++;
    }
  }
  return(M);
}



// [[Rcpp::export]]
NumericVector cal_vstar_h1_newsp(int n, int n_sample, sp_mat A, List nbs_list, IntegerVector dg){
  NumericVector M(n);
  int k=3;
  int src=0;
  while(src<n){
    if(dg(src)<1){
      M(src)=0;
      src++;
    }else{
      double s=0;
      IntegerVector nbs = nbs_list[src];
      int nid=nbs.size();
      IntegerVector visited(n);
      visited(src)=1;
      for(int ii=0; ii<nid; ii++){
        visited(nbs(ii))=1;
      }
      for(int j=1; j<n_sample; j++){
        arma::uvec S=get_sample_id_norep(nid,n-1);
        arma::uvec S_conc=get_sample_id_nov(nid,n-1,visited);
        arma::uvec indexes=sort(S);
        arma::uvec sorted=sort_index(S);
        int concurrent_ind=0;
        int i=0;
        while(i<nid){
          double t=0.0;
          int curr = indexes(i);
          int next=0;
          int nextnb=0;
          if(i+1 > (nid-1)){
             next = -1;
          }else{
             next = indexes(i+1);
             nextnb=nbs(sorted(i+1));
            }
          int currnb= nbs(sorted(i));
          
          if (next!=-1 && (curr%2==1 && next==curr+1)){
            if(A(currnb,nextnb)==0){
              t=1;}
              i=i+1;
          }else{
            int concurrent = S_conc(concurrent_ind);
            if (A(currnb,concurrent)==1){
              t=1;}
            concurrent_ind=concurrent_ind+1;
          }
          s=s+t;
          i=i+1;
        }
      }
      s=s/((n-1)/(k-1));
      M(src)=s/n_sample;
      src++;
    }
  }
  return(M);
}



// [[Rcpp::export]]
int in_same_bucket(IntegerVector SS){
  IntegerVector S=SS;
  int t=0;
  int b=3;
  if (SS.size()==3){
    if (S(0)%b==1&&S(1)==S(0)+1&&S(2)==S(1)+1){
    t=1;}
    else{
      t=0;}
    }
  if(SS.size()==2){
   if((S(0)%b==1&&S(1)==S(0)+1) or (S(0)%b==1&&S(1)==S(0)+2) or (S(0)%b==2&&S(1)==S(0)+1)){
   //if((S(0)%b==1&&S(1)==S(0)+1)){
    t=1;}
  else{
    t=0;}
}
  return(t);
  }




// [[Rcpp::export]]
NumericVector cal_3path_h1_newsp(int n, int n_sample, sp_mat A, List nbs_list, IntegerVector dg){
  NumericVector M(n);
  int k=4;
  int src=0;
  while(src<n){
    //Rcout << "src is" << src << std::endl;
    if(dg(src)<1){
      M(src)=0;
      src++;
    }else{
      double s=0;
      IntegerVector nbs = nbs_list[src];
      int nid=nbs.size();
      IntegerVector visited(n);
      visited(src)=1;
      for(int ii=0; ii<nid; ii++){
        visited(nbs(ii))=1;
      }
      for(int j=1; j<n_sample; j++){
        //Rcout << "j =" << j << std::endl;
        arma::uvec S=get_sample_id_norep(nid,n-1);
        int lenS_conc=2*nid;
        arma::uvec S_conc=get_sample_id_nov(lenS_conc,n-1,visited);
        //for(int ii=0; ii<lenS_conc; ii++){
        //  visited(S_conc(ii))=1;
        //}
        arma::uvec indexes=sort(S);
        arma::uvec sorted=sort_index(S);
        int concurrent_ind=0;
        int i=0;
        while(i<nid){
          double t=0.0;
          int curr = indexes(i);
          int next=0;
          int nextnext=0;
          int nextnb=0;
          int nextnextnb=0;
          if (i+1>(nid-1)){
            next=-1;
            nextnext=-1;}
          else{
            next = indexes(i+1);
            nextnb=nbs(sorted(i+1));
            if(i+2>(nid-1)){
              nextnext=-1;}
            else{
              nextnext=indexes(i+2);
              nextnextnb=nbs(sorted(i+2));}
          }
          int currnb= nbs(sorted(i));
            //int K2=dg(currnb);
            //int triangle_exist=0;
            IntegerVector SS1={curr,next,nextnext};
            IntegerVector SS2={curr,next};
            if(nextnext!=-1&& in_same_bucket(SS1)==1){
               t=0;
               i=i+2;}
            else if(next!=-1&&in_same_bucket(SS2)==1){
              int concurrent = S_conc(concurrent_ind);
              //S1=[src currnb nextnb concurrent];
              //t=check_path(A(S1,S1));
              t=(1-A(currnb,concurrent))*A(nextnb,concurrent)*(1-A(currnb,nextnb)) + (1-A(nextnb,concurrent))*A(currnb,concurrent)*(1-A(currnb,nextnb));
              i=i+1;
              concurrent_ind=concurrent_ind+1;}
            else{
              int concurrent = S_conc(concurrent_ind);
              int concurrent_nxt = S_conc(concurrent_ind+1);
            //S1=[src currnb concurrent concurrent_nxt];
            //t=check_path(A(S1,S1));
            t=A(currnb,concurrent)*A(concurrent,concurrent_nxt) *(1-A(currnb,concurrent_nxt)) + A(currnb,concurrent_nxt)*A(concurrent,concurrent_nxt)* (1-A(currnb,concurrent));
            concurrent_ind=concurrent_ind+2;
            }
            s=s+t;
            i=i+1;
        }
      }
      s=s/((n-1)/(k-1));
      M(src)=s/n_sample;
      src++;
    }
  }
  return(M);
}