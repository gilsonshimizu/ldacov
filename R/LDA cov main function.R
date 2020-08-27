#'Gibbs sampling LDA with covariates
#'
#'This function generates a sample of the parameters of interest using Slice sampling combined with Gibbs sampling.
#'
#'@param ncomm Integer. Number of clusters (communities, topics, etc.)
#'@param ngibbs Integer. Total sample size
#'@param nburn Integer. Beginning of the sample to be thrown away
#'@param y Matrix \eqn{L x S} with the quantities of each category in each sample unit. In text analysis we call this matrix bag of words.
#'@param xmat Matrix \eqn{L x D} with covariates for each sample unit.
#'@param phi.prior Numeric. Positive real number corresponding to the hyper parameter \eqn{\gamma} of the prior of \eqn{\phi}. We assume all \eqn{\gamma}'s equal.
#'@param array.lsk.init Array with initial values of the quantity of elements in each sample unit, category and cluster. Previously estimated in the model without covariates.
#'@param var.betas Vector. Positive real vector corresponding to the diagonal of the hyper parameter \eqn{T} of the prior of \eqn{\beta}.
#'@param phi.init Sample of the matrix \eqn{\Phi} previously generated with the model without covariates.
#'@param estimate.phi Logical. If FALSE \eqn{\Phi} will not be estimated. The matrix containing the posterior sample (\code{phi.init}) of the model without covariates will be used. If TRUE \eqn{\Phi} will be estimated. The matrix containing the posterior sample (\code{phi.init}) of the model without covariates will be used only to initialize the model.
#'@export
#'@return Returns a list containing samples of the parameters of interest (\code{phi}, \code{nlk} and \code{betas}) and a vector with the log likelihood (\code{llk}).
#'@examples
#'library(ldacov)
#'#Loads simulated data set
#'data("sim_data")
#'
#'#Model without covariates
#'res=gibbs.LDA(y=sim_data$y,
#'              ncomm=4,#Number of clusters found previously
#'              ngibbs=1000,
#'              nburn=500,
#'              psi=0.01,
#'              gamma=0.1)
#'
#'plot(res$llk,type='l')
#'
#'#Model with covariates
#'res.cov=gibbs.LDA.cov(ncomm=4,
#'                      ngibbs=1000,
#'                      nburn=500,
#'                      y=sim_data$y,
#'                      xmat=sim_data$xmat,
#'                      phi.prior=0.01,
#'                      array.lsk.init=res$array.lsk,
#'                      var.betas=c(10,rep(10,ncol(sim_data$xmat)-1)),
#'                      phi.init=res$phi,
#'                      estimate.phi=FALSE)
#'
#'plot(res.cov$llk,type='l')

gibbs.LDA.cov=function(ncomm,ngibbs,nburn,y,xmat,phi.prior,array.lsk.init,
                       var.betas,phi.init,estimate.phi){
  #basic settings
  nparam=ncol(xmat)
  nloc=nrow(y)
  nspp=ncol(y)
  ntot=apply(y,1,sum)

  #get phi by eliminating superfluous groups
  ncomm.init=ncol(phi.init)/nspp
  tmp=matrix(1:(ncomm.init*nspp),ncomm.init,nspp)
  seq1=1:ncomm
  ind1=tmp[-seq1,] #indicators for superfluous groups
  phi.mat=phi.init[,-ind1]
  phi.nrow=nrow(phi.mat)
  if (!estimate.phi) phi=matrix(phi.mat[1,],ncomm,nspp)
  if (estimate.phi)  phi=matrix(phi.mat[phi.nrow,],ncomm,nspp)

  #get theta
  nlk=apply(array.lsk.init,c(1,3),sum)
  theta1=nlk/apply(nlk,1,sum)

  #re-distribute individuals within array.lsk.init that are in eliminated communities
  seq1=1:ncomm
  array.lsk=array.lsk.init[,,seq1]
  for (i in 1:nloc){
    for (j in 1:nspp){
      tmp=array.lsk.init[i,j,-seq1]
      n=sum(tmp)
      if (n>0){
        prob=theta1[i,seq1]*phi[seq1,j]
        prob=prob/sum(prob)
        z=rmultinom(1,size=n,prob=prob)
        array.lsk[i,j,]=array.lsk[i,j,]+z
      }
    }
  }

  #initial values
  nlk=apply(array.lsk,c(1,3),sum)
  betas=matrix(0,nparam,ncomm)
  options(warn=-1) #sometimes I get "glm.fit: fitted rates numerically 0 occurred" here
  for (i in 1:ncomm){
    dat.tmp=cbind(nlk[,i],xmat[,-1])
    colnames(dat.tmp)=rep('',ncol(dat.tmp)) #this is important otherwise next line breaks when we only have a single covariate
    colnames(dat.tmp)[1]='y'
    dat.tmp1=as.data.frame(dat.tmp)
    res=try(glm.nb(y ~ ., data = dat.tmp1),silent=T)

    #if we run into an error using NB regression, use Poisson reg
    ind=grep('Error',res)
    if (length(ind)>0) res=glm(y~.,data=dat.tmp1,family='poisson')

    betas[,i]=res$coef
  }
  options(warn=2)
  nks=t(apply(array.lsk,2:3,sum))
  nk=rowSums(nks)
  # phi=nks/apply(nks,1,sum); apply(phi,1,sum)
  NBN=10

  #to store outcomes from gibbs sampler
  phi.out=matrix(NA,ngibbs,nspp*ncomm)
  nlk.out=matrix(NA,ngibbs,nloc*ncomm)
  llk.out=rep(NA,ngibbs)
  fmodel.out=matrix(NA,ngibbs,1)
  betas.out=matrix(NA,ngibbs,nparam*ncomm)
  NBN.out=matrix(NA,ngibbs,1)

  #useful stuff for slice sampler algorithm
  w.betas=1
  w.NBN=10
  MaxIter=100 #to avoid overly long slice samplers

  #to avoid numerical issues when calculating log(p) or log(1-p)
  LoThresh=0.00000001
  UpThresh=1-LoThresh

  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)

    #sample NBN
    media=exp(xmat%*%betas) #get mean
    NBN=SampleNBN(Media=media,y=nlk,NBN=NBN,w=w.NBN,MaxIter=MaxIter,LoThresh=LoThresh)

    #sample betas
    betas=SampleBetas(param=betas,y=nlk,xmat=xmat,w=w.betas,nparam=nparam,
                      ncomm=ncomm,var1=var.betas,NBN=NBN,MaxIter=MaxIter,
                      LoThresh=LoThresh)

    #sample phi
    if (estimate.phi)  phi=rdirichlet1(alpha=nks+phi.prior,ncomm=ncomm,nspp=nspp)
    if (!estimate.phi){
      oo=sample(phi.nrow,size=1)
      phi=matrix(phi.mat[oo,],ncomm,nspp)
    }
    # phi=phi.true

    #sample z
    media=exp(xmat%*%betas) #get mean
    NBP=NBN/(media+NBN)     #get NBP
    NBP[NBP>UpThresh]=UpThresh
    tmp = SampleArray(Arraylsk=array.lsk, nloc=nloc,nspp=nspp,ncomm=ncomm,NBN=NBN,
                      y=y,LogPhi=log(phi),LogOneMinusP=log(1-NBP),
                      runif1=runif(sum(y)),nlk=nlk)
    array.lsk=tmp$ArrayLSK
    # array.lsk=array.lsk.true

    nlk=apply(array.lsk,c(1,3),sum)
    # nks=t(apply(array.lsk,2:3,sum))
    nk=rowSums(nks)

    #calculate NB probabilities
    media.tmp=media
    media.tmp[media.tmp<LoThresh]=LoThresh
    p1=sum(dnbinom(nlk,mu=media.tmp,size=NBN,log=T))

    #calculate Multinom probabilities
    phi.tmp=phi
    phi.tmp[phi.tmp<LoThresh]=LoThresh
    tmp=LogLikMultin(nloc=nloc,ncomm=ncomm,nspp=nspp,LogPhi=log(phi.tmp), Arraylsk=array.lsk)
    p2=sum(tmp)

    #get phi prior
    p3=ldirichlet(x=phi.tmp,alpha=phi.prior)
    # log(ddirichlet(phi.tmp[2,],rep(phi.prior,nspp)))

    #get betas prior
    var.betas1=matrix(var.betas,nparam,ncomm)
    p4=dnorm(betas,mean=0,sd=sqrt(var.betas1),log=T)

    #store results
    llk.out[i]=sum(p1)+sum(p2)
    fmodel.out[i]=sum(p1)+sum(p2)+sum(p3)+sum(p4)
    phi.out[i,]=phi
    nlk.out[i,]=nlk
    betas.out[i,]=betas
    NBN.out[i]=NBN
  }

  list(llk=llk.out,phi=phi.out,nlk=nlk.out,betas=betas.out,fmodel=fmodel.out,NBN=NBN.out)
}


