#this function generates vmat, which is then used to generate the theta matrix
get.theta=function(nlk,gamma,ncomm,nloc,i,nburn,nks,theta,change1){
  #re-order groups from largest to smallest
  if (i%%change1 == 0 & i<nburn){
    med=colMeans(theta)
    ind=order(med,decreasing=T)
    nlk=nlk[,ind]
    nks=nks[ind,]
  }
  
  ngreater1=ngreater(nlk,nloc,ncomm)
  tmp=rbeta(n=nloc*(ncomm-1),shape1=nlk[,-ncomm]+1,shape2=ngreater1[,-1]+gamma)
  vmat=cbind(matrix(tmp,nloc,ncomm-1),1)
  theta=convertVtoTheta(vmat,rep(1,nloc))
  list(vmat=vmat,theta=theta,nks=nks)
}
#this function generates dirichlet random variables (1 one for each row of alpha)
rdirichlet1=function(alpha,ncomm,nspp){
  tmp=matrix(rgamma(n=ncomm*nspp,alpha,1),ncomm,nspp)
  soma=matrix(rowSums(tmp),ncomm,nspp)
  tmp/soma
}
