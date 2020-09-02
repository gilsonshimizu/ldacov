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
double LogTargetBetas(NumericVector LogMediaMiss, NumericVector param, NumericVector y, NumericMatrix xmat, 
                 int target, NumericVector var1, double NBN, double LgammaNBN) {
  NumericVector LogMediaComp=LogMediaMiss+xmat(_,target)*param[target];
  NumericVector media=exp(LogMediaComp);
  NumericVector NBP=NBN/(media+NBN);

  NumericVector LogPrior=(-1/(2*var1))*(param*param); //assumes mean zero
  double res=sum(lgamma(y+NBN)+NBN*log(NBP)+y*log(1-NBP))+sum(LogPrior);
  return res;
}

//this function doubles the interval until we are outside the slice
// [[Rcpp::export]]
NumericVector DoublingBetas(NumericVector LogMediaMiss, double yslice, double w, NumericVector param,
                       NumericVector y, NumericMatrix xmat,int target, NumericVector var1,
                       double NBN, double LgammaNBN, int MaxIter){
  NumericVector ParamLo=clone(param);
  NumericVector ParamHi=clone(param);
  ParamLo[target]=ParamLo[target]-w*runif(1)[0];
  ParamHi[target]=ParamLo[target]+w;
  double ylo=LogTargetBetas(LogMediaMiss,ParamLo,y,xmat,target,var1,NBN,LgammaNBN);
  double yhi=LogTargetBetas(LogMediaMiss,ParamHi,y,xmat,target,var1,NBN,LgammaNBN);

  int oo=0;
  int problem=0;
  while((ylo>yslice) & (oo<MaxIter)){
    ParamLo[target]=ParamLo[target]-w;
    ylo=LogTargetBetas(LogMediaMiss,ParamLo,y,xmat,target,var1,NBN,LgammaNBN);
    oo=oo+1;
  }
  if (oo >= MaxIter){
    problem=1;
  }
  oo=0;
  while((yhi>yslice) & (oo<MaxIter)){
    ParamHi[target]=ParamHi[target]+w;
    yhi=LogTargetBetas(LogMediaMiss,ParamHi,y,xmat,target,var1,NBN,LgammaNBN);
    oo=oo+1;
  }
  if (oo >= MaxIter){
    problem=1;
  }
  
  NumericVector res(2);
  res[0]=ParamLo[target];
  res[1]=ParamHi[target];
  if (problem==1){
    res[0]=param[target];
    res[1]=param[target];
  }
  return res;
}

//this function shrinks the slice if samples are outside the slice. If sample is inside the slice, accept this sample
// [[Rcpp::export]]
double SampleEachParamBetas(NumericVector LogMediaMiss,NumericVector rango1,double yslice,NumericVector param,
                            NumericVector y,NumericMatrix xmat,int target,NumericVector var1,
                            double NBN,int MaxIter, double LoThresh) {
  NumericVector param_orig1=clone(param);
  double yfim=R_NegInf;
  double x=0;
  double DistLo;
  double DistHi;
  double diff1=rango1[1]-rango1[0];
  double LgammaNBN=lgamma(NBN);
  int oo=0;
  while ((yfim<yslice) & (diff1 > LoThresh) & (oo<MaxIter)){
    x=rango1[0]+diff1*runif(1)[0]; //sample uniformly within this range
    param[target]=x;
    yfim=LogTargetBetas(LogMediaMiss,param,y,xmat,target,var1,NBN,LgammaNBN);
    if (yfim<yslice){ //shrink the slice if x falls outside
      DistLo=abs(rango1[0]-x);
      DistHi=abs(rango1[1]-x);
      if (DistLo<DistHi) rango1[0]=x;
      if (DistLo>DistHi) rango1[1]=x;
      diff1=rango1[1]-rango1[0];
    }
    oo=oo+1;
  }
  if ((diff1 <= LoThresh) | (oo >= MaxIter)){
    x=param_orig1[target];
  }
  return x;
}

//this function performs matrix multiplication after excluding i-th column in xmat and i-th element in param
// [[Rcpp::export]]
NumericVector MatrixMultip(NumericMatrix xmat, NumericVector param, int nparam) {
  NumericVector res(xmat.nrow());
  for(int i=0;i<nparam;i++){
    res=res+xmat(_,i)*param[i];
  }
  return res;
}

//this function samples all parameters for all communities, using a slice sampler for each parameter
// [[Rcpp::export]]
NumericMatrix SampleBetas(NumericMatrix param, NumericMatrix y, NumericMatrix xmat, 
                          double w,int nparam,int ncomm, NumericVector var1,double NBN, int MaxIter,
                          double LoThresh) {
  NumericVector LogMediaMiss(xmat.nrow());
  NumericVector LogMedia(xmat.nrow());
  double upper1;
  double yslice;
  NumericVector rango1(2);
  double LgammaNBN=lgamma(NBN);
  
  for(int i=0;i<ncomm;i++){
    LogMedia=MatrixMultip(xmat,param(_,i),nparam);
    for (int j=0; j<nparam; j++){
      //calculate logmean without the effect of the j-th covariate
      LogMediaMiss=LogMedia-xmat(_,j)*param(j,i);
      
      //define upper bound
      upper1=LogTargetBetas(LogMediaMiss,param(_,i),y(_,i),xmat,j,var1,NBN,LgammaNBN);
      yslice=upper1-rexp(1)[0]; //method suggest by Neal 2003 to sample uniformly vertically
      
      //define slice  
      rango1=DoublingBetas(LogMediaMiss,yslice,w,param(_,i),y(_,i),xmat,j,var1,NBN,LgammaNBN,MaxIter); //find range by doubling window
      
      //sample this particular parameter
      param(j,i)=SampleEachParamBetas(LogMediaMiss,rango1,yslice,param(_,i),y(_,i),xmat,j,var1,NBN,MaxIter,
                 LoThresh); //sample within the defined range (rango1)
      
      //re-calculate logmean
      LogMedia=LogMediaMiss+xmat(_,j)*param(j,i);
    }
  }
  return param;
}

