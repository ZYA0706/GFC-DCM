# This function was used to estimate the model parameters when K=5 in simulation.
Est_FC_TDCM_5 <- function(simdata,s,conditions,rep){
  #' Args: 
  #' @param simdata simulated data
  #' @param s index of sth condition
  #' @param conditions a data frame includes all conditions
  #' @param rep index of repth repelicaiton
  cond <- conditions[s, ]
  # file path
  path <- paste0("est/est_", cond$sim_name, "_rep=", rep, ".csv")
  
  # obtain data from simdata------------------------------------
  Nperson <- simdata$Nperson;Nstatement <- simdata$Nstatement
  Npair <- simdata$Npair;Nblock <- simdata$Nblock
  BS <- simdata$BS;K <- simdata$K;C <- simdata$C
  R <- simdata$R;Q <- simdata$Q
  pair_in <- simdata$pair_in
  all.patterns <- simdata$all.patterns;Y <- simdata$Y
  True_alpha<-simdata$true_alpha;True_lamda <- simdata$lamda
  # MCMC------------------------------------------------------------------------
  
  ## mycode---------------------------------------------------------------------
  
  mycode <- nimble::nimbleCode({
    
    for(n in 1:Nperson){
      for(b in 1:Npair){
        for (k in 1:K){
          w1[n, b, k] <- pow(alpha[n, k],Q1[b, k])
          w2[n, b, k] <- pow(alpha[n, k],Q2[b, k])}
        # first statement in paired item
        s_eta1[n, b] <- lamda1[pair_in[b,1]] * w1[n, b, 1] + lamda2[pair_in[b,1]] * w1[n, b, 2] + lamda3[pair_in[b,1]] * w1[n, b, 3] + lamda4[pair_in[b,1]] * w1[n, b, 4] + lamda5[pair_in[b,1]] * w1[n, b, 5]
        s_eta2[n, b] <- lamda12[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 2] + lamda13[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 3] + lamda14[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 4] + lamda15[pair_in[b,1]] * w1[n, b, 1] * w1[n, b, 5] + lamda23[pair_in[b,1]] * w1[n, b, 2] * w1[n, b, 3] + lamda24[pair_in[b,1]] * w1[n, b, 2] * w1[n, b, 4] + lamda25[pair_in[b,1]] * w1[n, b, 2] * w1[n, b, 5]+ lamda34[pair_in[b,1]] * w1[n, b, 3] * w1[n, b, 4] + lamda35[pair_in[b,1]] * w1[n, b, 3] * w1[n, b, 5] + lamda45[pair_in[b,1]] * w1[n, b, 4] * w1[n, b, 5]
        
        # second statement in paired item
        k_eta1[n, b] <- lamda1[pair_in[b,2]] * w2[n, b, 1] + lamda2[pair_in[b,2]] * w2[n, b, 2] + lamda3[pair_in[b,2]] * w2[n, b, 3] + lamda4[pair_in[b,2]] * w2[n, b, 4] + lamda5[pair_in[b,2]] * w2[n, b, 5]
        k_eta2[n, b] <- lamda12[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 2] + lamda13[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 3]  + lamda14[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 4] + lamda15[pair_in[b,2]] * w2[n, b, 1] * w2[n, b, 5] + lamda23[pair_in[b,2]] * w2[n, b, 2] * w2[n, b, 3] + lamda24[pair_in[b,2]] * w2[n, b, 2] * w2[n, b, 4] + lamda25[pair_in[b,2]]* w2[n, b, 2] * w2[n, b, 5] + lamda34[pair_in[b,2]] * w2[n, b, 3] * w2[n, b, 4] + lamda35[pair_in[b,2]] * w2[n, b, 3] * w2[n, b, 5] + lamda45[pair_in[b,2]] * w2[n, b, 4] * w2[n, b, 5]
        
        
        logit(p[n, b]) <- ((lamda0[pair_in[b,1]] + s_eta1[n, b] + s_eta2[n, b]) - (lamda0[pair_in[b,2]] + k_eta1[n, b] + k_eta2[n, b]))
        
        Y[n, b] ~ dbern(p[n, b])}
      
      for(k in 1:K) {alpha[n, k] <- all.patterns[c[n], k]}
      
      c[n] ~ dcat(pai[1:C])}
    
    pai[1:C] ~ ddirch(delta[1:C])
    
    for(i in 1:Nstatement) {
      item_parameter[i, 1:16] ~ dmnorm(item_mu[1:16], item_den[1:16,1:16])
      lamda0[i] <- item_parameter[i, 1]
      xlamda1[i] <- item_parameter[i, 2]
      xlamda2[i] <- item_parameter[i, 3]
      xlamda3[i] <- item_parameter[i, 4]
      xlamda4[i] <- item_parameter[i, 5]
      xlamda5[i] <- item_parameter[i, 6]
      xlamda12[i] <- item_parameter[i, 7]
      xlamda13[i] <- item_parameter[i, 8]
      xlamda14[i] <- item_parameter[i, 9]
      xlamda15[i] <- item_parameter[i, 10]
      xlamda23[i] <- item_parameter[i, 11]
      xlamda24[i] <- item_parameter[i, 12]
      xlamda25[i] <- item_parameter[i, 13]
      xlamda34[i] <- item_parameter[i, 14]
      xlamda35[i] <- item_parameter[i, 15]
      xlamda45[i] <- item_parameter[i, 16]
      
      lamda1[i] <- xlamda1[i] * Q[i, 1]
      lamda2[i] <- xlamda2[i] * Q[i, 2]
      lamda3[i] <- xlamda3[i] * Q[i, 3]
      lamda4[i] <- xlamda4[i] * Q[i, 4]
      lamda5[i] <- xlamda5[i] * Q[i, 5]
      lamda12[i]<-xlamda12[i] * Q[i,1] * Q[i,2]
      lamda13[i]<-xlamda13[i] * Q[i,1] * Q[i,3]
      lamda14[i]<-xlamda14[i] * Q[i,1] * Q[i,4]
      lamda15[i]<-xlamda15[i] * Q[i,1] * Q[i,5]
      lamda23[i]<-xlamda23[i] * Q[i,2] * Q[i,3]
      lamda24[i]<-xlamda24[i] * Q[i,2] * Q[i,4]
      lamda25[i]<-xlamda25[i] * Q[i,2] * Q[i,5]
      lamda34[i]<-xlamda34[i] * Q[i,3] * Q[i,4]
      lamda35[i]<-xlamda35[i] * Q[i,3] * Q[i,5]
      lamda45[i]<-xlamda45[i] * Q[i,4] * Q[i,5]
    }
    
    item_mu[1]~ dnorm(-2, 2)
    item_mu[2] ~ dnorm(4, 2)
    item_mu[3] ~ dnorm(4, 2)
    item_mu[4] ~ dnorm(4, 2)
    item_mu[5] ~ dnorm(4, 2)
    item_mu[6] ~ dnorm(4, 2)
    item_mu[7] ~ dnorm(0, 2)
    item_mu[8] ~ dnorm(0, 2)
    item_mu[9] ~ dnorm(0, 2)
    item_mu[10] ~ dnorm(0, 2)
    item_mu[11] ~ dnorm(0, 2)
    item_mu[12] ~ dnorm(0, 2)
    item_mu[13] ~ dnorm(0, 2)
    item_mu[14] ~ dnorm(0, 2)
    item_mu[15] ~ dnorm(0, 2)
    item_mu[16] ~ dnorm(0, 2)
    item_den[1:16,1:16] ~ dwish(R[1:16,1:16], 16)
    # Sigma_item[1:16,1:16] <- inverse(item_den[1:16,1:16])
  })
  
  ## myconstants----------------------------------------------------------------
  
  myconstants <- list(Nperson = Nperson,
                      Npair = Npair,
                      Nstatement = Nstatement,
                      Q = Q,
                      Q1 = Q[pair_in[,1],],
                      Q2 = Q[pair_in[,2],],
                      K = K,
                      C = C,
                      pair_in = pair_in,
                      delta = rep(1,C),
                      R = R)
  
  
  ## mydata---------------------------------------------------------------------
  
  mydata <- list(Y = Y,all.patterns = all.patterns)
  
  ## myinits--------------------------------------------------------------------
  
  gene_inits <- function(Nperson,Nstatement,C,R){
    item_mu <- c(-2,rep(4,log(C,2)),rep(0,choose(log(C,2),2)))
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
  ## mymonitors-----------------------------------------------------------------
  
  mymonitors <- c("alpha",
                  "lamda0","lamda1","lamda2","lamda3","lamda4","lamda5",
                  "lamda12","lamda13","lamda14","lamda15","lamda23",
                  "lamda24","lamda25","lamda34","lamda35","lamda45"
  )
  
  
  # output <- nimbleMCMC(code = mycode,constants = myconstants,data = mydata,
  #                      inits = myinits,monitors = mymonitors,thin = 1,
  #                      niter = 10000,nburnin = 5000,nchains = 1)
  
  ## run mcmc: parallerization--------------------------------------------------
  
  library(nimble)
  library(parallel)
  nbcores <- 2#detectCores() - 1
  my_cluster <- makeCluster(nbcores)
  
  ### mcmc workflow
  workflow <- function(seed, mydata, mycode,myinit,myconstants,mymonitors,niter,nburnin,thin) {
    
    library(nimble)
    
    mycode <- mycode
    
    # set.seed(123) # for reproducibility
    
    survival <- nimbleModel(code = mycode,
                            data = mydata,
                            inits = myinit,
                            constants = myconstants)
    Csurvival <- compileNimble(survival)
    Consurvival <- configureMCMC(survival,
                                 monitors = mymonitors,
                                 thin = thin)
    survivalMCMC <- buildMCMC(Consurvival)
    CsurvivalMCMC <- compileNimble(survivalMCMC,project = Csurvival)
    
    samples <- runMCMC(mcmc = CsurvivalMCMC,
                       niter = niter,
                       nburnin = nburnin,
                       setSeed = seed,thin = thin)
    
    return(samples)
  }
  
  X = sample(1:10000,2,F)
  
  output <- parLapply(cl = my_cluster,
                      X = X,# set the seed in workflow
                      fun = workflow,
                      mydata = mydata,
                      mycode = mycode,
                      myinit = gene_inits(Nperson = Nperson,Nstatement = Nstatement,C,R = R),
                      myconstants = myconstants,
                      mymonitors = mymonitors,
                      niter = 15000,nburnin = 7500,thin = 1)
  stopCluster(my_cluster)
  
  ## output---------------------------------------------------------------------
  library(MCMCvis)
  samples <- MCMCsummary(output)
  ### save posterior samples
  write.csv(samples, file = path, row.names = TRUE,col.names = TRUE)
  ### obtain a summary of the estimates all parameters
  alphaCols <- grep("alpha", rownames(samples))
  lamda0Cols <- grep("lamda0", rownames(samples))
  lamda1Cols <- grep("lamda1", rownames(samples))[1:Nstatement]
  lamda2Cols <- grep("lamda2", rownames(samples))[1:Nstatement]
  lamda3Cols <- grep("lamda3", rownames(samples))[1:Nstatement]
  lamda4Cols <- grep("lamda4", rownames(samples))[1:Nstatement]
  lamda5Cols <- grep("lamda5", rownames(samples))[1:Nstatement]
  lamda12Cols <- grep("lamda12", rownames(samples))
  lamda13Cols <- grep("lamda13", rownames(samples))
  lamda14Cols <- grep("lamda14", rownames(samples))
  lamda15Cols <- grep("lamda15", rownames(samples))
  lamda23Cols <- grep("lamda23", rownames(samples))
  lamda24Cols <- grep("lamda24", rownames(samples))
  lamda25Cols <- grep("lamda25", rownames(samples))
  lamda34Cols <- grep("lamda34", rownames(samples))
  lamda35Cols <- grep("lamda35", rownames(samples))
  lamda45Cols <- grep("lamda45", rownames(samples))
  
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
  est_lamda5 <- samples[lamda5Cols,1]
  est_lamda12 <- samples[lamda12Cols,1]
  est_lamda13 <- samples[lamda13Cols,1]
  est_lamda14 <- samples[lamda14Cols,1]
  est_lamda15 <- samples[lamda15Cols,1]
  est_lamda23 <- samples[lamda23Cols,1]
  est_lamda24 <- samples[lamda24Cols,1]
  est_lamda25 <- samples[lamda25Cols,1]
  est_lamda34 <- samples[lamda34Cols,1]
  est_lamda35 <- samples[lamda35Cols,1]
  est_lamda45 <- samples[lamda45Cols,1]
  est_lamda <- cbind(est_lamda0,
                     est_lamda1,est_lamda2,est_lamda3,est_lamda4,est_lamda5,
                     est_lamda12,est_lamda13,est_lamda14,est_lamda15,est_lamda23,
                     est_lamda24,est_lamda25,est_lamda34,est_lamda35,est_lamda45)
  
  
  ### parameter recovery
  
  PCCR <- sum(rowSums(est_alpha==True_alpha)==K)/Nperson
  ACCR <-  (colSums(est_alpha==True_alpha))/Nperson
  
  Esti_Intercept <- para_recovery(est_lamda0,True_lamda[,1])
  
  Esti_Main <- para_recovery(est_lamda[,2:6],True_lamda[,2:6])

  Esti_Interaction <- para_recovery(est_lamda[,7:16],True_lamda[,7:16])
  # Save results----------------------------------------------------------------
  data <- new.env()
  data$PCCR <- PCCR
  data$ACCR <- ACCR
  data$Esti_Intercept <- Esti_Intercept
  data$Esti_Main <- Esti_Main
  data$Esti_Interaction <- Esti_Interaction
  data <- as.list(data)
  return(data)
}
