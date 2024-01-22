################################################################################
#                                                                              #
#     Generalized Forced-Choice Diagnostic Classification Model Based on       #
#                                                                              #
#                   Thurstone's Law of Comparative Judgement                   #
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


# Step:2 MCMC------------------------------------------------------------------------

## mycode-----------------------------------------------------------------------

mycode <- nimble::nimbleCode({
  
  for(n in 1:Nperson){
    for(b in 1:Npair){
      for (k in 1:K){
        w1[n, b, k] <-pow(alpha[n, k],Q1[b, k])
        w2[n, b, k] <- pow(alpha[n, k],Q2[b, k])}
      # first statement in paired item
      s_eta1[n, b] <- lamda1[pair_in[b,1]] * w1[n, b, 1] + lamda2[pair_in[b,1]] * w1[n, b, 2] + lamda3[pair_in[b,1]] * w1[n, b, 3] + lamda4[pair_in[b,1]] * w1[n, b, 4]
      s_eta2[n, b] <- lamda12[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 2] + lamda13[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 3] + lamda14[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 4] + lamda24[pair_in[b,1]] * w1[n, b, 2] * w1[n, b, 4] + lamda34[pair_in[b,1]] * w1[n, b, 3] * w1[n, b, 4]
      s_eta3[n, b] <- lamda234[pair_in[b,1]] * w1[n, b, 2] * w1[n, b, 3] * w1[n, b, 4]
      
      # second statement in paired item
      k_eta1[n, b] <- lamda1[pair_in[b,2]] * w2[n, b, 1] + lamda2[pair_in[b,2]] * w2[n, b, 2] + lamda3[pair_in[b,2]] * w2[n, b, 3] + lamda4[pair_in[b,2]] * w2[n, b, 4]
      k_eta2[n, b] <- lamda12[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 2] + lamda13[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 3]  + lamda14[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 4] + lamda24[pair_in[b,2]] * w2[n, b, 2] * w2[n, b, 4] + lamda34[pair_in[b,2]] * w2[n, b, 3] * w2[n, b, 4]
      k_eta3[n, b] <- lamda234[pair_in[b,2]] * w2[n, b, 2] * w2[n, b, 3] * w2[n, b, 4]
      
      
      logit(P[n,b]) <- ((lamda0[pair_in[b,1]] + s_eta1[n, b] + s_eta2[n, b] + s_eta3[n, b]) - (lamda0[pair_in[b,2]] + k_eta1[n, b] + k_eta2[n, b] + k_eta3[n, b]))
      
      Y[n, b] ~ dbern(P[n,b])}
    
    for(k in 1:K) {alpha[n, k] <- all.patterns[c[n], k]}
    
    c[n] ~ dcat(pai[1:C])}
  
  pai[1:C] ~ ddirch(delta[1:C])
  
  for(i in 1:Nstatement) {
    item_parameter[i, 1:11] ~ dmnorm(item_mu[1:11], item_den[1:11,1:11])
    lamda0[i] <- item_parameter[i, 1]
    xlamda1[i] <- item_parameter[i, 2]
    xlamda2[i] <- item_parameter[i, 3]
    xlamda3[i] <- item_parameter[i, 4]
    xlamda4[i] <- item_parameter[i, 5]
    xlamda12[i] <- item_parameter[i, 6]
    xlamda13[i] <- item_parameter[i, 7]
    xlamda14[i] <- item_parameter[i, 8]
    xlamda24[i] <- item_parameter[i, 9]
    xlamda34[i] <- item_parameter[i, 10]
    xlamda234[i] <- item_parameter[i, 11]
    
    lamda1[i] <- xlamda1[i]  *  Q[i, 1]
    lamda2[i] <- xlamda2[i]  *  Q[i, 2]
    lamda3[i] <- xlamda3[i]  *  Q[i, 3]
    lamda4[i] <- xlamda4[i]  *  Q[i, 4]
    lamda12[i] <- xlamda12[i]  *  Q[i,1]  *  Q[i,2]
    lamda13[i] <- xlamda13[i]  *  Q[i,1]  *  Q[i,3]
    lamda14[i] <- xlamda14[i]  *  Q[i,1]  *  Q[i,4]
    lamda24[i] <- xlamda24[i]  *  Q[i,2]  *  Q[i,4]
    lamda34[i] <- xlamda34[i]  *  Q[i,3]  *  Q[i,4]
    lamda234[i] <- xlamda24[i]  *  Q[i,2] * Q[i,3]  *  Q[i,4]
    
  }
  
  item_mu[1]~ dnorm(-2, 2)
  item_mu[2] ~ dnorm(3, 2)
  item_mu[3] ~ dnorm(3, 2)
  item_mu[4] ~ dnorm(3, 2)
  item_mu[5] ~ dnorm(3, 2)
  item_mu[6] ~ dnorm(0, 2)
  item_mu[7] ~ dnorm(0, 2)
  item_mu[8] ~ dnorm(0, 2)
  item_mu[9] ~ dnorm(0, 2)
  item_mu[10] ~ dnorm(0, 2)
  item_mu[11] ~ dnorm(0, 2)
  item_den[1:11,1:11] ~ dwish(R[1:11,1:11], 11)
  # Sigma_item[1:16,1:16] <- inverse(item_den[1:16,1:16])
  
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
                    Q = Q,
                    Q1 = Q[pair_in[,1],],
                    Q2 = Q[pair_in[,2],],
                    pair_in = pair_in,
                    K = K,
                    C = C,
                    delta = rep(1,C),
                    R = R)


## mydata------------------------------------------------------------------------

Y <- Y_binary; all.patterns <- GDINA::attributepattern(K)
mydata <- list(Y = Y,all.patterns = all.patterns)

## myinits-----------------------------------------------------------------------

gene_inits <- function(Nperson,Nstatement,C,R){
  item_mu <- c(-2,rep(3,log(C,2)),rep(0,choose(log(C,2),2)))
  item_den <- matrix(rWishart(1,df = nrow(R),Sigma = R),nrow = nrow(R),ncol = nrow(R))
  item_parameter <- matrix(0,nrow = Nstatement,ncol = nrow(R))
  for (i in 1:Nstatement) {
    item_parameter[i,] <- rmnorm_chol(n = 1,mean = item_mu,cholesky = item_den)
  }
  c <- sample(1:C,Nperson,replace = T)
  pai <- rep(1/C,C)
  return(list(item_mu = item_mu,
              item_den = item_den,
              item_parameter = item_parameter,
              c = c,
              pai = pai))
}

myinits <- gene_inits(Nperson = Nperson,Nstatement = Nstatement,C = C,R = R)

## mymonitors--------------------------------------------------------------------

mymonitors <- c("alpha",
                "lamda0","lamda1","lamda2","lamda3","lamda4",
                "lamda12","lamda13","lamda14",
                "lamda24","lamda34","lamda234",
                "ppp_test","ppp_pair")


## ## one-step MCMC----------------------------------------------------------------

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
WAIC_GFC_DCM <- output$WAIC

## parallerization-----------------------------

# nbcores <- detectCores() - 1
# my_cluster <- makeCluster(nbcores)

### mcmc workflow
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
#                     myinit = gene_inits(Nperson = Nperson,Nstatement = Nstatement,C,R = R),
#                     myconstants = myconstants,
#                     mymonitors = mymonitors,
#                     niter = 15000,nburnin = 7500,
#                     thin = 1)
# stopCluster(my_cluster)

# str(output)

## output------------------------------------------------------------------------

library(MCMCvis)
samples <- MCMCsummary(output$samples)
### obtain a summary of the estimates all parameters
alphaCols <- grep("alpha", rownames(samples))
lamda0Cols <- grep("lamda0", rownames(samples))
lamda1Cols <- grep("lamda1", rownames(samples))[1:Nstatement]
lamda2Cols <- grep("lamda2", rownames(samples))[1:Nstatement]
lamda3Cols <- grep("lamda3", rownames(samples))[1:Nstatement]
lamda4Cols <- grep("lamda4", rownames(samples))[1:Nstatement]
lamda12Cols <- grep("lamda12", rownames(samples))
lamda13Cols <- grep("lamda13", rownames(samples))
lamda14Cols <- grep("lamda14", rownames(samples))
lamda24Cols <- grep("lamda24", rownames(samples))
lamda34Cols <- grep("lamda34", rownames(samples))
lamda234Cols <- grep("lamda234", rownames(samples))
ppp_testCols <- grep("ppp_test", rownames(samples))
ppp_pairCols <- grep("ppp_pair", rownames(samples))

# View(est_lamda)

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

est_lamda0 <- samples[lamda0Cols,1]
est_lamda1 <- samples[lamda1Cols,1]
est_lamda2 <- samples[lamda2Cols,1]
est_lamda3 <- samples[lamda3Cols,1]
est_lamda4 <- samples[lamda4Cols,1]
est_lamda12 <- samples[lamda12Cols,1]
est_lamda13 <- samples[lamda13Cols,1]
est_lamda14 <- samples[lamda14Cols,1]
est_lamda24 <- samples[lamda24Cols,1]
est_lamda34 <- samples[lamda34Cols,1]
est_lamda234 <- samples[lamda234Cols,1]
est_lamda <- cbind(est_lamda0,
                   est_lamda1,est_lamda2,est_lamda3,est_lamda4,
                   est_lamda12,est_lamda13,est_lamda14,
                   est_lamda24,est_lamda34,est_lamda234)

est_ppp_test <- samples[ppp_testCols,1]
est_ppp_pair <- samples[ppp_pairCols,1]

psd_lamda0 <- round((samples[lamda0Cols,"sd"]),3)
psd_lamda1 <- round((samples[lamda1Cols,"sd"]),3)
psd_lamda2 <- round((samples[lamda2Cols,"sd"]),3)
psd_lamda3 <- round((samples[lamda3Cols,"sd"]),3)
psd_lamda4 <- round((samples[lamda4Cols,"sd"]),3)
psd_lamda12 <- round((samples[lamda12Cols,"sd"]),3)
psd_lamda13 <- round((samples[lamda13Cols,"sd"]),3)
psd_lamda14 <- round((samples[lamda14Cols,"sd"]),3)
psd_lamda24 <- round((samples[lamda24Cols,"sd"]),3)
psd_lamda34 <- round((samples[lamda34Cols,"sd"]),3)
psd_lamda234 <- round((samples[lamda234Cols,"sd"]),3)
psd_lamda <- c(psd_lamda0,psd_lamda1,psd_lamda2,psd_lamda3,psd_lamda4,
               psd_lamda12,psd_lamda13,psd_lamda14,psd_lamda24,psd_lamda34,
               psd_lamda234)

# Step 3: Save results----------------------------------------------------------------
GFC_DCM_LIST <- list(data = mydata,Q = Q,
                     Nperson = Nperson,Nblock = Nblock,
                     Npair = Npair,Nstatement = Nstatement,K = K,
                     samples = samples,
                     GFC_DCM_code = mycode,
                     WAIC = WAIC_GFC_DCM,
                     ppp_test = est_ppp_test,
                     ppp_pair = est_ppp_pair,
                     ppp_block = apply(matrix(est_ppp_pair,ncol = 3,byrow = T),1,mean),
                     est_alpha = est_alpha,
                     est_lamda = est_lamda,
                     psd_lamda = psd_lamda,
                     Rhat_all = samples[,c("Rhat")],
                     n.eff_all = samples[,c("n.eff")],
                     NO_CON_RATE = NO_CON_RATE
)

saveRDS(object = GFC_DCM_LIST,file = 'result/RD_GFC_DCM_list.rds')
saveRDS(object = output,file = 'result/RD_GFC_DCM_mcmcoutput.rds')
