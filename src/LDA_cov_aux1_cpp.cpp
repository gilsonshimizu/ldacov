// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
// [[Rcpp::export]]
int whichLessDVPresence(double value, NumericVector prob) {
  int res=prob.length()-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

// This function calculates the multinomial distribution
// [[Rcpp::export]]
double ldmultinom(NumericVector x, double size, NumericVector lprob) {
  NumericVector res1=x*lprob;
  double res=lgamma(size+1) - sum(lgamma(x+1)) + sum(res1);
  return(res);
}

// this function calculates the loglikelih based on multinomial
// [[Rcpp::export]]
NumericVector LogLikMultin(int nloc,int ncomm, int nspp, NumericMatrix LogPhi, 
                           NumericVector Arraylsk){
  NumericVector p2(nloc);
  //convert array into arma::cube
  NumericVector vecArray(Arraylsk);
  arma::cube ArrayLSK1(vecArray.begin(), nloc, nspp, ncomm, false);
  NumericVector tmp(nspp);
  double soma;
  for (int l = 0; l < nloc; l++) {
    for (int k = 0; k < ncomm; k++){
      for (int s = 0; s < nspp; s++){
        tmp[s]=ArrayLSK1(l,s,k);
      }
      soma=sum(tmp);
      if (soma>0) p2[l]=p2[l]+ldmultinom(tmp,soma,LogPhi(k,_));
    }
  }
  return(p2);
}

// This function samples Array.lsk
// [[Rcpp::export]]

List SampleArray(NumericVector Arraylsk, int nloc, int nspp, int ncomm, double NBN,
                 IntegerMatrix y, NumericMatrix LogPhi, NumericMatrix LogOneMinusP,
                 NumericVector runif1, IntegerMatrix nlk){
  
  //convert array into arma::cube
  NumericVector vecArray=clone(Arraylsk);
  arma::cube ArrayLSK1(vecArray.begin(), nloc, nspp, ncomm, false);
  arma::cube ArrayLSK1Orig(Arraylsk.begin(), nloc, nspp, ncomm, false);

  //initialize stuff
  NumericVector prob(ncomm);
  int ind;
  int oo=0;
  
  for (int l = 0; l < nloc; l++) {
    //go over each individual (i.e., each element of ArrayLSK1Orig)
    for (int s = 0; s < nspp; s++){
      if (y(l,s)>0){
        for (int k = 0; k < ncomm; k++){
          if (ArrayLSK1Orig(l,s,k)>0){
            for (int i = 0; i < ArrayLSK1Orig(l,s,k); i++){
              //remove i-th individual
              ArrayLSK1(l,s,k)=ArrayLSK1(l,s,k)-1;
              nlk(l,k)=nlk(l,k)-1;

              //calculate assignment probabilities
              for (int k1 = 0; k1 < ncomm; k1++){
                prob[k1]=log(nlk(l,k1)+NBN)-log(ArrayLSK1(l,s,k1)+1)+LogPhi(k1,s)+LogOneMinusP(l,k1);  
              }
              prob=prob-max(prob);
              prob=exp(prob);
              prob=prob/sum(prob);
              
              //sample cluster membership of i-th individual
              ind=whichLessDVPresence(runif1[oo], prob);
              oo=oo+1;
              
              //update counts
              ArrayLSK1(l,s,ind)=ArrayLSK1(l,s,ind)+1;
              nlk(l,ind)=nlk(l,ind)+1;
            }
          }
        }
      }
    }
  }

  List L = List::create(Named("ArrayLSK") =ArrayLSK1);
  
  return(L);
}
