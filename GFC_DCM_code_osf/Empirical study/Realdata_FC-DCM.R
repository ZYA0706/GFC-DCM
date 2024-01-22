################################################################################
#                                                                              #
#            Forced-Choice Diagnostic Classification Model                     #
#                                                                              #
#                           Empirical Study                                    #
#                                                                              #
################################################################################

# Notes：
# Please create the corresponding folder first to ensure that the program can run normally.
# Empirical_study>>realdata
# Empirical_study>>result

# Step:1 Preparation------------------------------------------------------------
working='GFC-DCM_code_osf/Empirical study'
setwd(working)
rm(list = ls())
gc()
## packages
library(GDINA)
library(gtools)
library(nimble)
library(parallel)

## load real dataset
Y_rank <- read.csv(file = 'realdata/FC-MPDQ.csv',header = F,sep = ',')
for (i in 1:nrow(Y_rank)) {
  for (j in 1:ncol(Y_rank)) {
    if (Y_rank[i,j] == 4) {
      Y_rank[i,j] <- 1
    }else if(Y_rank[i,j] == 3){
      Y_rank[i,j] <- 2
    }else if(Y_rank[i,j] == 2){
      Y_rank[i,j] <- 3
    }else if(Y_rank[i,j] == 1){
      Y_rank[i,j] <- 4
    }
  }
}

Nperson <- nrow(Y_rank)           # NO.person
Nstatement <- ncol(Y_rank)        # NO.statements
BS <- 4                           # NO.statements per block
Nblock <- Nstatement / BS         # NO.block
Npair <- Nblock  *  choose(BS,2)  # NO.pair
K <- 4                            # NO.attribute
C <- 2^K                          # NO.class
R <- diag(11)


### Index for statements in all blocks
ori_block <- matrix(data = c(3,33,36,40,8,9,11,13,5,6,37,42,12,19,28,41,20,22,27,35,
                             15,16,18,24,25,30,43,44,1,10,31,38,7,14,21,29,17,26,32,34,
                             2,4,23,39),nrow = Nblock,byrow = T)
### Q-matrix
ori_Q <- read.csv(file = 'realdata/MPDQ-Q.csv',header = F,sep = ',')[1:44,]
Q <- ori_Q[c(t(ori_block)),]

### Index for all paired items

block <- matrix(1:Nstatement,nrow = Nblock,byrow = T)
pair_in <- NULL
for (i in 1:Nblock) {
  pair_in <- append(x = pair_in,as.vector(combn(x = block[i,],2)))
}
pair_in <- matrix(pair_in,ncol = 2,byrow = T)

### binary response data of all paired items for all respondents
Y_binary <- matrix(0,nrow = Nperson,ncol = Nblock * choose(BS,2)) 
for (n in 1:Nperson) {
  i_r <- matrix(0,nrow = Nblock,ncol = choose(BS,2))
  for (m in 1:Nblock) {
    temp_combm <- t(combn(block[m,],2)) 
    #------------------------
    for (t in 1:nrow(temp_combm)) {
      judge <- Y_rank[n,temp_combm[t,]]
      
      if (judge[1] > judge[2]) {
        i_r[m,t] <- 1
      }else{
        i_r[m,t] <- 0
      }
    }
  }
  Y_binary[n,] <- as.vector(t(i_r))
}
# MCMC------------------------------------------------------------------------

## mycode-----------------------------------------------------------------------


# data

mycode <- nimble::nimbleCode({
  
  for(n in 1:Nperson){
    for(b in 1:Npair){
      for (k in 1:K){
        w1[n, b, k] <- pow(alpha[n, k],Q1[b, k])
        w2[n, b, k] <- pow(alpha[n, k],Q2[b, k])}
      
      P[n,b] <-  d0[b] + (0.5 - d0[b]) * step(prod(w1[n, b, 1:K]) - prod(w2[n, b, 1:K])) + (sd[Npair_in[b,1]] * prod(w1[n, b, 1:K]) + sd[Npair_in[b,2]] * (1 - prod(w2[n, b, 1:K]))) * step(prod(w1[n, b, 1:K]) - prod(w2[n, b, 1:K]) - 1)
      # P[n,b] <-  d0[b] + (d1[b] - d0[b]) * step(prod(w1[n, b, 1:K]) - prod(w2[n, b, 1:K])) + (sd[Npair_in[b,1]] * prod(w1[n, b, 1:K]) + sd[Npair_in[b,2]] * (1 - prod(w2[n, b, 1:K]))) * step(prod(w1[n, b, 1:K]) - prod(w2[n, b, 1:K]) - 1)
      
      Y[n,b] ~ dbern(P[n,b])
      
    }
  }
  
  for (n in 1:Nperson) {
    for (k in 1:K) {
      logit(att_prob[n, k])<-1.7*delta1[k]*(theta[n]-delta0[k])
      alpha[n, k] ~ dbern(att_prob[n, k])}
    theta[n] ~ dnorm(0,1)
  }
  
  for (k in 1:K) {
    delta1[k] ~ dlnorm(0, 1)
    delta0[k] ~ dnorm(0, 0.25)
  }
  
  for (b in 1:Npair) {
    d0[b] ~ T(dbeta(1,1), ,0.5)
    # d1[b] ~ T(dbeta(1,1), 0.4, 0.6)
  }
  
  for(i in 1:Nstatement) {
    sd[i] ~ T(dbeta(1,1), ,0.5)
  }
  
  
  
  # posterior predictive p-value:ppp
  for (n in 1:Nperson) {
    for (b in 1:Npair) {
      
      # predictive value
      y_rep[n,b] ~ dbern(P[n,b])
      
      # calculate standardized residual
      
      res_rep[n,b] <- pow(y_rep[n,b] - P[n,b],2) / (P[n,b]  *  (1 - P[n,b]))
      
      res_data[n,b] <- pow(Y[n,b] - P[n,b],2) / (P[n,b]  *  (1 - P[n,b]))
      
      
    }
  }
  
  # test ppp
  
  sum_res_rep <- sum(res_rep[1:Nperson,1:Npair])
  sum_res_data <- sum(res_data[1:Nperson,1:Npair])
  ppp_test <- nimStep(sum_res_rep - sum_res_data)
  
  # pair ppp
  
  for (b in 1:Npair) {
    sum_pair_rep[b] <- sum(res_rep[1:Nperson,b])
    sum_pair_data[b] <- sum(res_data[1:Nperson,b])
    ppp_pair[b] <- nimStep(sum_pair_rep[b] - sum_pair_data[b])
  }
  
})


## myconstants-------------------------------------------------------------------

myconstants <- list(Nperson = Nperson,
                    Npair = Npair,
                    Nstatement = Nstatement,
                    # Q = Q,
                    Q1 = Q[Npair_in[,1],],
                    Q2 = Q[Npair_in[,2],],
                    Npair_in = Npair_in,
                    K = K)


## mydata------------------------------------------------------------------------

Y <- Y_binary
mydata <- list(Y = Y)

## myinits-----------------------------------------------------------------------

gene_inits <- function(Nperson,Nstatement,Npair,K){
  d0 <- runif(Npair,0,0.5)
  # d1 <- runif(Npair,0.4,0.6)
  sd <- runif(Nstatement,0,0.5)
  theta <- rep(0,Nperson)
  delta1 <- rlnorm(K,0,1)
  delta0 <- rnorm(K,0,0.25)
  return(list(d0 = d0,
              # d1 = d1,
              sd = sd,
              theta = theta,
              delta1 = delta1,
              delta0 = delta0
  ))
}

myinits <- gene_inits(Nperson = Nperson,Npair = Npair,Nstatement = Nstatement,K = K)

## mymonitors--------------------------------------------------------------------

mymonitors <- c("alpha","theta",
                "d0","sd",
                # "d1",
                "delta1","delta0",
                "ppp_test","ppp_pair")

# one-step ---------------------------------------------------------
start<-proc.time()

output <- nimble::nimbleMCMC(code = mycode,
                             constants = myconstants,
                             data = mydata,
                             inits = myinits,
                             monitors = mymonitors,
                             niter = 20000,
                             nburnin = 10000,
                             nchains = 2,
                             thin = 1,
                             WAIC = T)

time<-(proc.time()-start)
WAIC_FC_DCM <- output$WAIC

## parallerization-----------------------------

# library(nimble)
# library(parallel)
# nbcores <- detectCores() - 1
# my_cluster <- makeCluster(nbcores)
# 
# mcmc workflow
# workflow <- function(seed, mydata, mycode,myinit,myconstants,mymonitors,niter,nburnin,thin) {
#   
#   library(nimble)
#   
#   mycode <- mycode
#   
#   # set.seed(123) # for reproducibility
#   
#   survival <- nimbleModel(code = mycode,
#                           data = mydata,
#                           inits = myinit,
#                           constants = myconstants)
#   Csurvival <- compileNimble(survival)
#   Consurvival <- configureMCMC(survival,
#                                monitors = mymonitors,
#                                thin = thin)
#   survivalMCMC <- buildMCMC(Consurvival)
#   CsurvivalMCMC <- compileNimble(survivalMCMC,project = Csurvival)
#   
#   samples <- runMCMC(mcmc = CsurvivalMCMC,
#                      niter = niter,
#                      nburnin = nburnin,
#                      setSeed = seed,thin = thin)
#   
#   return(samples)
# }
# X <- sample(1:10000,2)
# output <- parLapply(cl = my_cluster,
#                     X = X,# set the seed in workflow
#                     fun = workflow,
#                     mydata = mydata,
#                     mycode = mycode,
#                     myinit = gene_inits(Nperson = Nperson,Nstatement = Nstatement,Npair = Npair,K = K),
#                     myconstants = myconstants,
#                     mymonitors = mymonitors,
#                     niter = 15000,nburnin = 7500,
#                     thin = 1)
# stopCluster(my_cluster)
# 
# # str(output)




## output------------------------------------------------------------------------

library(MCMCvis)
samples <- MCMCsummary(output$samples)
### obtain a summary of the estimates all parameters
alphaCols <- grep("alpha", rownames(samples))
thetaCols <- grep("theta", rownames(samples))
d0Cols <- grep("d0", rownames(samples))
d1Cols <- grep("d1", rownames(samples))[1:Nstatement]
sdCols <- grep("sd", rownames(samples))[1:Nstatement]
delta1Cols <- grep("delta1", rownames(samples))
delta0Cols <- grep("delta0", rownames(samples))
ppp_testCols <- grep("ppp_test", rownames(samples))
ppp_pairCols <- grep("ppp_pair", rownames(samples))

Rhat_convergence <- samples[,c("Rhat")]
no_convergence <- rownames(samples)[which(Rhat_convergence > 1.1)]
NO_CON_RATE <- length(no_convergence) / length(Rhat_convergence)
# View(samples[,c("Rhat","n.eff")])

est_alpha <- matrix(data = samples[alphaCols,1],nrow = Nperson,K)
for (i in 1:nrow(est_alpha)) {
  for (j in 1:K) {
    if (est_alpha[i,j] >= 0.5) {
      est_alpha[i,j] <- 1
    }else(
      est_alpha[i,j] <- 0
    )
  }
}


est_theta <- samples[thetaCols,1]
psd_theta <- round(samples[thetaCols,"sd"],3)


est_d0 <- samples[d0Cols,1]
est_d1 <- samples[d1Cols,1]
est_sd <- samples[sdCols,1]
psd_d0 <- round(samples[d0Cols,"sd"],3)
psd_d1 <- round(samples[d1Cols,"sd"],3)
psd_sd <- round(samples[sdCols,"sd"],3)
est_ppp_test <- samples[ppp_testCols,1]
est_ppp_pair <- samples[ppp_pairCols,1]


est_delta1 <- samples[delta1Cols,1]
est_delta0 <- samples[delta0Cols,1]
psd_delta1 <- round(samples[delta1Cols,"sd"],3)
psd_delta0 <- round(samples[delta0Cols,"sd"],3)

# Step 3: Save results----------------------------------------------------------
FC_DCM_LIST <- list(data = mydata,Q = Q,
                    Nperson = Nperson,Nblock = Nblock,
                    Npair = Npair,Nstatement = Nstatement,K = K,
                    samples = samples,
                    FC_TDCM_code = mycode,
                    WAIC = WAIC_FC_DCM,
                    ppp_test = est_ppp_test,
                    ppp_pair = est_ppp_pair,
                    ppp_block = apply(matrix(est_ppp_pair,ncol = 3,byrow = T),1,mean),
                    est_P = est_P, 
                    est_alpha = est_alpha,
                    est_theta = est_theta,
                    est_d0 = est_d0,
                    est_d1 = est_d1,
                    est_sd = est_sd,
                    est_delta0 = est_delta0,
                    est_delta1 = est_delta1,
                    psd_theta = psd_theta,
                    psd_d0 = psd_d0,
                    psd_d1 = psd_d1,
                    psd_sd = psd_sd,
                    psd_delta0 = psd_delta0,
                    psd_delta1 = psd_delta1,
                    Rhat_all = samples[,c("Rhat")],
                    n.eff_all = samples[,c("n.eff")],
                    NO_CON_RATE = NO_CON_RATE
)

saveRDS(object = FC_DCM_LIST,file = 'result/RD_FC_DCM_list.rds')
saveRDS(object = output,file = 'result/RD_FC_DCM_mcmcoutput.rds')


