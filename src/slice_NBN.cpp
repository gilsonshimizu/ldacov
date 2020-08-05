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

// This function calculates the logtarget distribution
// [[Rcpp::export]]
double LogTargetNBN(NumericMatrix Media, NumericMatrix y, double NBN) {
  NumericVector NBP;
  double res=0;
  NumericVector tmp;
  for (int j = 0; j < Media.ncol(); j++){
    NBP=NBN/(NBN+Media(_,j));  
    tmp=lgamma(y(_,j)+NBN)-lgamma(NBN)-lgamma(y(_,j)+1)+NBN*log(NBP)+y(_,j)*log(1-NBP);
    res=res+sum(tmp);
  }
  return res;
}

//this function doubles the interval until we are outside the slice
// [[Rcpp::export]]
NumericVector DoublingNBN(double yslice, double w, NumericMatrix y, double NBN,NumericMatrix Media,int MaxIter,
                          double LoThresh){
  double ParamLo=NBN;
  double ParamHi=NBN;
  ParamLo=ParamLo-w*runif(1)[0];
  if (ParamLo<LoThresh) ParamLo=LoThresh;
  ParamHi=ParamLo+w;
  double ylo=LogTargetNBN(Media,y,ParamLo);
  double yhi=LogTargetNBN(Media,y,ParamHi);

  int oo=0;
  while((ylo>yslice) & (ParamLo>0) & (oo<MaxIter)){
    ParamLo=ParamLo-w;
    if (ParamLo<LoThresh) ParamLo=LoThresh;
    ylo=LogTargetNBN(Media,y,ParamLo);
    oo=oo+1;
  }
  oo=0;
  while((yhi>yslice) & (ParamHi<1000) & (oo<MaxIter)){
    ParamHi=ParamHi+w;
    yhi=LogTargetNBN(Media,y,ParamHi);
    oo=oo+1;
  }
  NumericVector res(2);
  res[0]=ParamLo;
  res[1]=ParamHi;
  return res;
}

//this function shrinks the slice if samples are outside the slice. If sample is inside the slice, accept this sample
// [[Rcpp::export]]
double ShrinkAndSample(NumericMatrix Media,NumericVector rango1,double yslice,NumericMatrix y,double NBN, int MaxIter,
                       double LoThresh) {
  double yfim=R_NegInf;
  double x=NBN;
  double DistLo;
  double DistHi;
  double diff1=rango1[1]-rango1[0];
  int oo=0;
  while ((yfim<yslice) & (diff1 > LoThresh) & (oo<MaxIter)){
    x=rango1[0]+diff1*runif(1)[0]; //sample uniformly within this range
    yfim=LogTargetNBN(Media,y,x);
    if (yfim<yslice){ //shrink the slice if x falls outside
      DistLo=abs(rango1[0]-x);
      DistHi=abs(rango1[1]-x);
      if (DistLo<DistHi) rango1[0]=x;
      if (DistLo>DistHi) rango1[1]=x;
      diff1=rango1[1]-rango1[0];
    }
    oo=oo+1;
  }
  return x;
}

//this function samples NBN using a slice sampler 
// [[Rcpp::export]]
double SampleNBN(NumericMatrix Media,NumericMatrix y,double NBN, double w, int MaxIter, double LoThresh) {
  double upper1;
  double yslice;
  NumericVector rango1(2);

  //define upper bound
  upper1=LogTargetNBN(Media,y,NBN);
  yslice=upper1-rexp(1)[0]; //method suggest by Neal 2003 to sample uniformly vertically
      
  //define slice  
  rango1=DoublingNBN(yslice, w, y, NBN, Media, MaxIter,LoThresh); //find range by doubling window
  
  //sample this particular parameter
  NBN=ShrinkAndSample(Media,rango1,yslice,y,NBN, MaxIter,LoThresh) ; //sample within the defined range (rango1)
      
  return NBN;
}
