#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector randomShuffle(IntegerVector a) {

  std::random_shuffle(a.begin(), a.end(), randWrapper);
  
  return a;
}


// [[Rcpp::export]]
NumericMatrix cal_tri_h1(int n, int n_sample, NumericMatrix A){
  NumericMatrix M(n,n_sample);
  double t1=0.0;
  int k=3;
  for(int src=0; src<n; src++){
    t1=0;
    IntegerVector ni= seq(0,n-1);
    ni.erase(src);
    for(int j=0; j<n_sample; j++){
      IntegerVector id = randomShuffle(ni);
      double t=0;
      for(int ii=0; ii<floor((n-1)/(k-1)); ii++){
        IntegerVector S=id[seq((ii*(k-1)),((ii+1)*(k-1)-1))];
        t=t+A(src,S[0])*A(S[0],S[1])*A(S[1],src);
    }
      t=t/((n-1)/(k-1));
      t1=t1+t;
      M(src,j)=t1/j;
  }
  }
  return(M);
}
/* If want to return one number Tnhat should just do return(M[,N]) the last column */


// [[Rcpp::export]]
NumericMatrix cal_vstar_h1(int n, int n_sample, NumericMatrix A){
  NumericMatrix M(n,n_sample);
  double t1=0.0;
  int k=3;
  for(int src=0; src<n; src++){
    t1=0;
    IntegerVector ni= seq(0,n-1);
    ni.erase(src);
    for(int j=0; j<n_sample; j++){
      IntegerVector id = randomShuffle(ni);
      double t=0;
      for(int ii=0; ii<floor((n-1)/(k-1)); ii++){
        IntegerVector S=id[seq((ii*(k-1)),((ii+1)*(k-1)-1))];
        t=t+A(S[1],S[0])*A(S[0],src)*(1-A(S[1],src)) + A(S[0],S[1])*A(S[1],src)*(1-A(S[0],src)) + A(S[0],src)*A(src,S[1])*(1-A(S[0],S[1]));
      }
      t=t/((n-1)/(k-1));
      t1=t1+t;
      M(src,j)=t1/j;
    }
  }
  return(M);
}

// [[Rcpp::export]]
NumericMatrix cal_4cycle_h1(int n, int n_sample, NumericMatrix A){
  NumericMatrix M(n,n_sample);
  double t1=0.0;
  int k=4;
  for(int src=0; src<n; src++){
    t1=0;
    IntegerVector ni= seq(0,n-1);
    ni.erase(src);
    for(int j=0; j<n_sample; j++){
      IntegerVector id = randomShuffle(ni);
      double t=0;
      for(int ii=0; ii<floor((n-1)/(k-1)); ii++){
        IntegerVector S=id[seq((ii*(k-1)),((ii+1)*(k-1)-1))];
        t=t+(A(S[0],S[1])*A(S[1],S[2])*A(S[2],src)*A(src,S[0])*(1-A(S[0],S[2]))*(1-A(S[1],src))
               + A(S[0],S[1])*A(S[1],src)*A(src,S[2])*A(S[2],S[0])*(1-A(S[0],src))*(1-A(S[1],S[2]))
               + A(S[0],S[2])*A(S[2],S[1])*A(S[1],src)*A(src,S[0])*(1-A(S[0],S[1]))*(1-A(S[2],src)) );
      }
      t=t/((n-1)/(k-1));
      t1=t1+t;
      M(src,j)=t1/j;
    }
  }
  return(M);
}
/* If want to return one number Tnhat should just do return(M[,N]) the last column */


