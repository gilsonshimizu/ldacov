# ldacov
ldacov is an R package to estimate an LDA model with covariates using MCMC algorithms.
## Installation
Use the devtools package to install ldacov directly from github.
```
install_github("gilsonshimizu/ldacov")
```
## Usage
```
library(ldacov)
#Loads simulated data set
data("sim_data")

#Model without covariates
res=gibbs.LDA(y=sim_data$y,
             ncomm=4,#Number of clusters found previously
             ngibbs=1000,
             nburn=500,
             psi=0.01,
             gamma=0.1)

plot(res$llk,type='l')

#Model with covariates
res.cov=gibbs.LDA.cov(ncomm=4,
                     ngibbs=1000,
                     nburn=500,
                     y=sim_data$y,
                     xmat=sim_data$xmat,
                     phi.prior=0.01,
                     array.lsk.init=res$array.lsk,
                     var.betas=c(10,rep(10,ncol(sim_data$xmat)-1)),
                     phi.init=res$phi)

plot(res.cov$llk,type='l')
```
