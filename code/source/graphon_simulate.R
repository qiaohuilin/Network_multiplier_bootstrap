library(methods)
library(Rcpp)
library(RcppArmadillo)

# Define  objects
setClass("graphon",representation(w_function="character"))
setClass("sbm",representation(B_mat="matrix",pi_vec="numeric"))
setClass("lsm",representation(kernel_function="character",latent_dist ="character"))

#Some notes on proper syntax for the classes above:
# For graphon objects, the w_function should be a string that contains a Rcpp function, which takes as
# input 2 double variables and returns a double variable.  At this point, sanity checks are not
# included.  Note that the function must be named "w_function".  

# For lsm objects, the kernel_function is defined in an identical manner to the w_function for graphons. The latent_dist slot
# should be a character vector consisting of two components.  These components correspond to the characters
# used before and after the number of observations to be simulated (n) is specified, in Rcpp syntax.
#For example, if the latent distribution is N(0,1), then latent_dist = c('rnorm(' , ',0,1)') is appropriate.
# For now, the latent positions must be double variables.  

setClass("sparse_graphon",contains = "graphon", representation  = list(nu_seq = "function"),prototype(nu_seq = function(x) return(1)))
setClass("sparse_sbm",contains = "sbm",representation  = list(nu_seq = "function"),prototype(nu_seq = function(x) return(1)))
setClass("sparse_lsm",contains = "lsm",representation  = list(nu_seq = "function"),prototype(nu_seq = function(x) return(1)))

sparse_object = function(x,...) UseMethod("sparse_object")

sparse_object.sbm = function(sbm,nu_seq) {
  output = new("sparse_sbm", B_mat =sbm@B_mat, pi_vec =sbm@pi_vec, nu_seq= nu_seq)
  return(output)
}

sparse_object.graphon = function(graphon,nu_seq) {
  output = new("sparse_graphon", w_function =graphon@w_function, nu_seq= nu_seq)
  return(output)
}


# Functions below generate a function that take as input a sparse object and produces a Rcpp function to simulate a graph of size n from the model.  For now, the Rcpp code is included inline; I will consider moving the code to a .cpp file.

create_simulation_function =  function(x) UseMethod("create_simulation_function")

create_simulation_function.sparse_sbm = function(sparse_sbm){
  
  cppFunction( 'arma::mat my_function( int n, double nu, NumericMatrix B_mat,
               NumericVector pi_vec ) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               int num_class = pi_vec.size();
               IntegerVector my_classes = seq_len(num_class);
               
               IntegerVector class_assignment = Rcpp::sample(my_classes,n,true,pi_vec);
               
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj < n; jj++) 
               {
               int x_temp = class_assignment[ii]-1;
               int y_temp = class_assignment[jj]-1;
               
               double P_temp = nu * B_mat(x_temp, y_temp);
               A_mat(ii,jj) = R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
               }
               ',depends = "RcppArmadillo")
  
  
  output = function(n) {
    return( my_function(n=n,nu=sparse_sbm@nu_seq(n),B_mat = sparse_sbm@B_mat, pi_vec = sparse_sbm@pi_vec))
  }
  return(output)
}

create_simulation_function.sparse_graphon = function(sparse_graphon) {
  
  cppFunction( 'arma::mat my_function( int n, double nu) {
                            arma::mat A_mat(n,n,arma::fill::zeros);
                            NumericVector latent_positions = runif(n,0,1);
                            for (int ii = 0; ii < n; ii++)
                            {
                            for(int jj = ii+1; jj< n; jj++) 
                            {
                            double P_temp = nu * w_function(latent_positions[ii],latent_positions[jj]);
                            A_mat(ii,jj) =  R::rbinom(1,P_temp);
                            }
                            }
                            A_mat = arma::symmatu(A_mat);
                            return(A_mat);
}
',depends = "RcppArmadillo", includes = sparse_graphon@w_function)

  
  output = function (n) {
    return(my_function(n=n, sparse_graphon@nu_seq(n)))
  }
  return(output)
}

############################

# Graphon Function Repository, selected from Chan and Airoldi (2014)


#Comments: eigenvalues/functions of graphon calculable following Xu (2017), high rank.  
graphon.1 = 'double w_function( double u, double v) {
  return(std::abs(u-v));
}'

my.graphon.1 = new("graphon", w_function = graphon.1)

graphon.2 = 'double w_function( double u, double v) {
  return((pow(u,2)+pow(v,2))/3*cos(1/(pow(u,2)+pow(v,2)))+0.15);
}'

my.graphon.2 = new("graphon", w_function = graphon.2)


###########################################
# SBM Repository

#my.sbm = new("sbm",B_mat= matrix(c(0.4,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.7),nrow=3), pi_vec=c(0.3,0.3,0.4))

my.sbm = new("sbm",B_mat= matrix(c(0.6,0.2,0.2,0.2),nrow=2), pi_vec=c(0.65,0.35))

my.sbm1 = new("sbm",B_mat= matrix(c(0.6,0.2,0.2,0.2),nrow=2), pi_vec=c(0.5,0.5))
