#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
IntegerVector findAdjacentStates2(sp_mat adjacency, int col) {
  IntegerVector out;
  sp_mat::const_col_iterator start = adjacency.begin_col(col);
  sp_mat::const_col_iterator end = adjacency.end_col(col);
  for ( sp_mat::const_col_iterator i = start; i != end; ++i )
  {
    out.push_back(i.row());
  }
  return out;
}

// [[Rcpp::export]]
List find_all(int n,sp_mat A){
  List l(n);
  for(int i=0;i<n;i++){
   IntegerVector nbs = findAdjacentStates2(A,i);
   l[i]=nbs;
  }
  return l ;
}

