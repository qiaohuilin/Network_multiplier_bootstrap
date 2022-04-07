#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int count_4cycles(int n, NumericMatrix A){
  int s=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      for(int k=0; k<j; k++){
        for(int t=0; t<k; t++){
          s=s+(A(i,j)*A(j,k)*A(k,t)*A(t,i)*(1-A(i,k))*(1-A(j,t))
                 + A(i,j)*A(j,t)*A(t,k)*A(k,i)*(1-A(i,t))*(1-A(j,k))
                 + A(i,k)*A(k,j)*A(j,t)*A(t,i)*(1-A(i,j))*(1-A(k,t)));
             }
      }
    }
  }
  return s;
}

// [[Rcpp::export]]
NumericVector count_4cycles_h1(int n, NumericMatrix A){
  NumericVector s(n);
  for(int i=0; i<n; i++){
    IntegerVector ni= seq(0,n-1);
    ni.erase(i);
    for(int jj=0; jj<n-1; jj++){
      for(int kk=0; kk<jj; kk++){
        for(int tt=0; tt<kk; tt++){
          int j=ni[jj];
          int k=ni[kk];
          int t=ni[tt];
          s[i]=s[i]+(A(i,j)*A(j,k)*A(k,t)*A(t,i)*(1-A(i,k))*(1-A(j,t))
                       + A(i,j)*A(j,t)*A(t,k)*A(k,i)*(1-A(i,t))*(1-A(j,k))
                       + A(i,k)*A(k,j)*A(j,t)*A(t,i)*(1-A(i,j))*(1-A(k,t)));
        }
      }
    }
  }
  return s;
}


// [[Rcpp::export]]
NumericMatrix count_4cycles_h2(int n, NumericMatrix A){
  NumericMatrix s(n,n);
  for(int i=0; i<n; i++){
    //IntegerVector ni= seq(0,n-1);
    //ni.erase(i);
    for(int j=0; j<i;j++){
      //ni.erase(j);
      for(int k=0; k<n; k++){
        if(k!=i & k!=j){
        for(int t=0; t<k; t++){
          if(t!=i & t!=j){
          //int k=kk;
          //int t=tt;
          //int k=ni[kk];
          //int t=ni[tt];
          s(i,j)=s(i,j)+(A(i,j)*A(j,k)*A(k,t)*A(t,i)*(1-A(i,k))*(1-A(j,t))
                           + A(i,j)*A(j,t)*A(t,k)*A(k,i)*(1-A(i,t))*(1-A(j,k))
                           + A(i,k)*A(k,j)*A(j,t)*A(t,i)*(1-A(i,j))*(1-A(k,t)));
          }}
        }}
      s(j,i)=s(i,j);
    }
  }
  return s;
}

// [[Rcpp::export]]
int count_4cliques(int n, NumericMatrix A){
  int s=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      for(int k=0; k<j; k++){
        for(int t=0; t<k; t++){
          s=s+A(i,j)*A(j,k)*A(k,t)*A(t,i)*A(i,k)*A(j,t);
        }
      }
    }
  }
  return s;
}


// [[Rcpp::export]]
int count_3path(int n, NumericMatrix A){
  int s=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(j != i){
        for(int k=0; k<n; k++){
          if(k != i and k!=j){
            for(int t=0; t<n; t++){
              if(t !=i and t!=j and t!=k){
                s=s+A(i,j)*A(j,k)*A(k,t)*(1-A(i,t))*(1-A(j,t))*(1-A(i,k));
              }}
          }}
      }}
  }
  s=s/2;
  return (s);
}


// [[Rcpp::export]]
NumericVector count_3path_i(int n, NumericMatrix A){
  NumericVector s(n);
  for(int i=0; i<n; i++){
    int ss=0;
    for(int j=0; j<n; j++){
      if(j != i){
        for(int k=0; k<n; k++){
          if(k != i and k!=j){
            for(int t=0; t<n; t++){
              if(t !=i and t!=j and t!=k){
                ss=ss+A(i,j)*A(j,k)*A(k,t)*(1-A(i,t))*(1-A(j,t))*(1-A(i,k))+ A(j,i)*A(i,k)*A(k,t)*(1-A(j,t))*(1-A(i,t))*(1-A(j,k));
              }}
          }}
      }}
    s(i)=ss;
  }
  return (s);
}
