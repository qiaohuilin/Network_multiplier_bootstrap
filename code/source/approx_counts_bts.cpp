#include <random>
#include <vector>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp ;


/* sample without replacement  */

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
NumericVector cal_kato_tri_h1(int n, int n_sample, NumericMatrix A){
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
      sum= sum + (A(i,id3[0])*A(id3[0],id3[1])*A(id3[1],i));
    }
    g1[i]=sum;
  }
  return(g1);
}




// [[Rcpp::export]]
NumericVector cal_kato_vstar_h1(int n, int n_sample, NumericMatrix A){
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
      sum= sum + (A(i,id3[0])*A(id3[0],id3[1])*(1-A(id3[1],i))) + (A(i,id3[1])*A(id3[1],id3[0])*(1-A(id3[0],i)))+(A(id3[0],i)*A(i,id3[1])*(1-A(id3[0],id3[1]))); ;
    }
    g1[i]=sum;
  }
  return(g1);
}

// [[Rcpp::export]]
NumericVector cal_kato_4cycle_h1(int n, int n_sample, NumericMatrix A){
  NumericVector g1(n);
  for(int i=0; i<n; i++){
    IntegerVector nn= seq(0,n-1);
    nn.erase(i);
    double sum=0.0;
    for(int iter=0; iter<n_sample; iter++){
      int k=3;
      arma::uvec id = get_sample_id_norep(k,n);
      IntegerVector id3(k);
      for (int ii=0; ii<k; ii++) {
        int jj = id(ii) -1 ;  // arma 
        id3[ii] = nn[jj]; // templated
      }
      sum= sum +  (A(id3[0],id3[1])*A(id3[1],id3[2])*A(id3[2],i)*A(i,id3[0])*(1-A(id3[0],id3[2]))*(1-A(id3[1],i))
                     + A(id3[0],id3[1])*A(id3[1],i)*A(i,id3[2])*A(id3[2],id3[0])*(1-A(id3[0],i))*(1-A(id3[1],id3[2]))
                     + A(id3[0],id3[2])*A(id3[2],id3[1])*A(id3[1],i)*A(i,id3[0])*(1-A(id3[0],id3[1]))*(1-A(id3[2],i)) );
    }
    g1[i]=sum;
  }
  return(g1);
}

