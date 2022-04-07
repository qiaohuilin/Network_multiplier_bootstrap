#include <random>
#include <vector>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp ;

// [[Rcpp::export]]
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
IntegerVector randomShuffle(IntegerVector a) {
  
  std::random_shuffle(a.begin(), a.end(), randWrapper);
  
  return a;
}


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
NumericVector calcpdgg1kato(int n, int n_sample, NumericMatrix A, NumericVector label){
  NumericVector g1(n);
  for(int i=0; i<n; i++){
    IntegerVector nn= seq(0,n-1);
    nn.erase(i);
    double sum=0.0;
    for(int iter=0; iter<n_sample; iter++){
      int k=1;
      arma::uvec id = get_sample_id_norep(k,n);
      IntegerVector id3(k);
      for (int ii=0; ii<k; ii++) {
        int jj = id(ii) -1 ;  // arma 
        id3[ii] = nn[jj]; // templated
      }
      if(label[i] != label[id3[0]]){
        sum= sum + A(i,id3[0]);
      }
    }
    g1[i]=sum;
  }
  return(g1);
}

// [[Rcpp::export]]
NumericVector calcptrig1kato(int n, int n_sample, NumericMatrix A, NumericVector label){
  NumericVector g1(n);
  for(int i=0; i<n; i++){
    IntegerVector nn= seq(0,n-1);
    nn.erase(i);
    double sum=0.0;
    for(int iter=0; iter<n_sample; iter++){
      int k=2;
      arma::uvec id = get_sample_id_norep(k,n);
      IntegerVector id3(k);
      for (int ii=0; ii<k; ii++) {
        int jj = id(ii) -1 ;  // arma 
        id3[ii] = nn[jj]; // templated
      }
      if(label[i]+label[id3[0]]+label[id3[1]] != 300 & label[i]+label[id3[0]]+label[id3[1]] != 600){
        sum= sum + (A(i,id3[0])*A(id3[0],id3[1])*A(id3[1],i));
      }
    }
    g1[i]=sum;
  }
  return(g1);
}

// [[Rcpp::export]]
NumericMatrix calcpdgg1(int n, int n_sample, NumericMatrix A, NumericVector label){
  NumericMatrix M(n,n_sample);
  double t1=0.0;
  int k=2;
  for(int src=0; src<n; src++){
    t1=0;
    IntegerVector ni= seq(0,n-1);
    ni.erase(src);
    for(int j=0; j<n_sample; j++){
      IntegerVector id = randomShuffle(ni);
      double t=0;
      for(int ii=0; ii<floor((n-1)/(k-1)); ii++){
        IntegerVector S=id[seq((ii*(k-1)),((ii+1)*(k-1)-1))];
        if(label[src] != label[S[0]]){
          t=t+A(src,S[0]);
        }
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
NumericMatrix calcptrig1(int n, int n_sample, NumericMatrix A, NumericVector label){
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
        if(label[src]+label[S[0]]+label[S[1]] != 300 & label[src]+label[S[0]]+label[S[1]] != 600){
          t=t+A(src,S[0])*A(S[0],S[1])*A(S[1],src);
        }
      }
      t=t/((n-1)/(k-1));
      t1=t1+t;
      M(src,j)=t1/j;
    }
  }
  return(M);
}
/* If want to return one number Tnhat should just do return(M[,N]) the last column */
