#'Gibbs sampling LDA without covariates
#'
#'This function generates a sample of the parameters of interest using Metropolis Hastings algorithm combined with Gibbs sampling.
#'
#'@param y Matrix \eqn{L x S} with the quantities of each category in each sample unit. In text analysis we call this matrix bag of words.
#'@param ncomm Integer. Number of clusters (communities, topics, etc.)
#'@param ngibbs Integer. Total sample size
#'@param nburn Integer. Beginning of the sample to be thrown away
#'@param psi Numeric. Positive real number corresponding to the hyper parameter of the prior of \eqn{\phi}.
#'@param gamma Numeric. Positive real number corresponding to the hyper parameter of the prior of \eqn{V}.
#'@export

LDA.abundance=function(y,ncomm,ngibbs,nburn,psi,gamma){
  #get data
  nspp=ncol(y)
  nloc=nrow(y)

  #useful stuff
  hi=0.999999
  lo=0.000001

  #initial values of parameters
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(1/nspp,ncomm,nspp)

  #gibbs details
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  vmat.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  llk=rep(NA,ngibbs)
  # log.prior=rep(NA,ngibbs)

  options(warn=2)
  zeroes=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:ngibbs){
    print(i)

    #re-order z from time to time
    if (i<nburn & i%%50==0){
      med=apply(theta,2,median)
      ordem=order(med,decreasing=T)
      theta=theta[,ordem]
      phi=phi[ordem,]
    }

    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp, zeroes=zeroes)
    array.lsk=tmp$ArrayLSK
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true

    #get parameters
    tmp=get.theta(nlk=nlk,gamma=gamma,ncomm=ncomm,nloc=nloc,i=i,nburn=nburn,nks=nks,theta=theta,change1=20)
    nks=tmp$nks
    vmat=tmp$vmat
    theta=tmp$theta
    # theta[theta>hi]=hi; theta[theta<lo]=lo
    # theta=theta.true

    phi=rdirichlet1(alpha=nks+psi,ncomm=ncomm,nspp=nspp)
    # phi[phi>hi]=hi; phi[phi<lo]=lo
    # phi=phi.true

    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo

    #calculate log prior (I often get Inf!!)
    # vmat1=vmat[,-ncomm]
    # vmat1[vmat1>hi]=hi; vmat1[vmat1<lo]=lo
    # log.p.betas=sum(dbeta(vmat1,1,gamma,log=T))
    # log.p.phi=sum(log(ddirichlet(phi,rep(psi,nspp))))
    # print(c(log.p.betas,log.p.phi))

    #store results
    llk[i]=sum(y*log(prob))
    # log.prior[i]=log.p.betas+log.p.phi
    theta.out[i,]=theta
    phi.out[i,]=phi
    vmat.out[i,]=vmat
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       # log.prior=log.prior[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,],
       vmat=vmat.out[seq1,],
       array.lsk=array.lsk)
}
