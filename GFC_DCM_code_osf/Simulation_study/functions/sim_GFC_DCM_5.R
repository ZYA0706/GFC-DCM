# This function was used to generate the data when K=5 in simulation.
sim_FC_TDCM_5 <- function(s,conditions,rep){
  #' Args: 
  #' @param s index of sth condition
  #' @param conditions a data frame includes all conditions
  #' @param rep index of repth repelicaiton
  # simulation condition--------------------------------------------------------
  
  cond <- conditions[s, ]
  # file path
  path <- paste0("simdata/sim_", cond$sim_name, "_rep=", rep, ".txt")
  
  cat(cond$sim_name, "_rep=", rep, sep = "", "\n")
  
  Nperson <- cond$Npersons              # NO.person
  Nstatement <- cond$Nstatement         # NO.statements
  BS <- cond$BS                         # NO.statements per block
  K <- cond$K                           # NO.attributes
  LSM <- cond$LSM                       # latent strutrue model
  # random seed
  set.seed(Nperson * 100 + Nstatement + rep)
  
  Nblock <- Nstatement / BS             # NO.block
  Npair <- Nblock * choose(BS,2)        # NO.pair
  C <- 2^K                              # NO.class
  R <- diag(16)
  block <- develop_block2345(item_number = 1:Nstatement,Block_size = BS)
  
  # simulation data-------------------------------------------------------------
  ## Q-matrix-------------------------------------------------------------------
  library(GDINA)
  all.patterns <- attributepattern(K)
  att_1 <- all.patterns[which(rowSums(all.patterns) == 1),]
  att_2 <- all.patterns[which(rowSums(all.patterns) == 2),]
  Q <- matrix(NA, Nstatement, K)
  Q[1:(Nstatement/2), ] <- att_1[sample(1:nrow(att_1), (Nstatement/2), replace = T), ]
  Q[(Nstatement/2 + 1 ):(Nstatement), ] <- att_2[sample(1:nrow(att_2), (Nstatement/2), replace = T), ]
  Q[1:(Nstatement),] <- Q[sample(1:nrow(Q),Nstatement,replace = F),]
  colnames(Q) <- colnames(all.patterns)
  Q <- Q[c(t(block)),]
  
  new_block <- matrix(1:Nstatement,Nblock,BS,byrow = T)
  pair_in <- matrix(ncol = 2)
  for (b in 1:Nblock) {
    pair_in <- rbind(pair_in,t(combn(x = new_block[b,],2)))
  }
  pair_in <- pair_in[-1,]
  
  ## person parameters-------------------------------------- -------------------
  
  if (LSM == 'unif') {
    ### uniform distribution
    true_alpha <- matrix(NA,Nperson,K)
    for (n in 1:Nperson) {
      true_alpha[n,] <- all.patterns[sample(1:C,1,T),]
    }
  }else if(LSM == 'mvnorm'){
    ### multivariate normal distribution
    gs <- data.frame(guess=c(runif(Nstatement, 0, 0.2)),slip=c(runif(Nstatement, 0, 0.2)))
    m <- rep(0, K)
    # cov.alpha <- matrix(0, K, K)
    cov.alpha <- matrix(0.5, K, K)
    diag(cov.alpha) <- 1
    cutoffs <- qnorm(c(1: K) / (K + 1))
    sim_GDINA <- simGDINA(Nperson, Q, gs.parm = gs,
                          att.dist = "mvnorm",gs.args=list("mono.constraint"),
                          mvnorm.parm = list(mean = m, sigma = cov.alpha, cutoffs = cutoffs),
                          model = "GDINA")
    true_alpha <- sim_GDINA$attribute
  }
  
  
  
  
  
  
  ## statement parameters-------------------------------------------------------
  
  ###  the means of statement parameters
  
  xlamda0_mu <- -2
  xlamda1_mu <- 4
  xlamda2_mu <- 4
  xlamda3_mu <- 4
  xlamda4_mu <- 4
  xlamda5_mu <- 4
  xlamda12_mu <- 0
  xlamda13_mu <- 0
  xlamda14_mu <- 0
  xlamda15_mu <- 0
  xlamda23_mu <- 0
  xlamda24_mu <- 0
  xlamda25_mu <- 0
  xlamda34_mu <- 0
  xlamda35_mu <- 0
  xlamda45_mu <- 0
  
  item_mu <- c(xlamda0_mu,xlamda1_mu,xlamda2_mu,xlamda3_mu,xlamda4_mu,xlamda5_mu,
               xlamda12_mu,xlamda13_mu,xlamda14_mu,xlamda15_mu,xlamda23_mu,xlamda24_mu,
               xlamda25_mu,xlamda34_mu,xlamda35_mu,xlamda45_mu)
  
  ### covariance matrix
  
  item_sigma <- diag(16)
  diag(item_sigma) <- 0.5
  item_mvrnorm <- MASS::mvrnorm(Nstatement, mu = item_mu, Sigma = item_sigma, 16)
  xlamda0 <- as.matrix(item_mvrnorm[,1])
  xlamda1 <- as.matrix(item_mvrnorm[,2])
  xlamda2 <- as.matrix(item_mvrnorm[,3])
  xlamda3 <- as.matrix(item_mvrnorm[,4])
  xlamda4 <- as.matrix(item_mvrnorm[,5])
  xlamda5 <- as.matrix(item_mvrnorm[,6])
  xlamda12 <- as.matrix(item_mvrnorm[,7])
  xlamda13 <- as.matrix(item_mvrnorm[,8])
  xlamda14 <- as.matrix(item_mvrnorm[,9])
  xlamda15 <- as.matrix(item_mvrnorm[,10])
  xlamda23 <- as.matrix(item_mvrnorm[,11])
  xlamda24 <- as.matrix(item_mvrnorm[,12])
  xlamda25 <- as.matrix(item_mvrnorm[,13])
  xlamda34 <- as.matrix(item_mvrnorm[,14])
  xlamda35 <- as.matrix(item_mvrnorm[,15])
  xlamda45 <- as.matrix(item_mvrnorm[,16])
  
  lamda1<-matrix(NA,Nstatement,1)
  lamda2<-matrix(NA,Nstatement,1)
  lamda3<-matrix(NA,Nstatement,1)
  lamda4<-matrix(NA,Nstatement,1)
  lamda5<-matrix(NA,Nstatement,1)
  lamda12<-matrix(NA,Nstatement,1)
  lamda13<-matrix(NA,Nstatement,1)
  lamda14<-matrix(NA,Nstatement,1)
  lamda15<-matrix(NA,Nstatement,1)
  lamda23<-matrix(NA,Nstatement,1)
  lamda24<-matrix(NA,Nstatement,1)
  lamda25<-matrix(NA,Nstatement,1)
  lamda34<-matrix(NA,Nstatement,1)
  lamda35<-matrix(NA,Nstatement,1)
  lamda45<-matrix(NA,Nstatement,1)
  
  for(i in 1:Nstatement){
    lamda1[i]<-xlamda1[i]*Q[i,1]
    lamda2[i]<-xlamda2[i]*Q[i,2]
    lamda3[i]<-xlamda3[i]*Q[i,3]
    lamda4[i]<-xlamda4[i]*Q[i,4]
    lamda5[i]<-xlamda5[i]*Q[i,5]
    lamda12[i]<-xlamda12[i]*Q[i,1]*Q[i,2]
    lamda13[i]<-xlamda13[i]*Q[i,1]*Q[i,3]
    lamda14[i]<-xlamda14[i]*Q[i,1]*Q[i,4]
    lamda15[i]<-xlamda15[i]*Q[i,1]*Q[i,5]
    lamda23[i]<-xlamda23[i]*Q[i,2]*Q[i,3]
    lamda24[i]<-xlamda24[i]*Q[i,2]*Q[i,4]
    lamda25[i]<-xlamda25[i]*Q[i,2]*Q[i,5]
    lamda34[i]<-xlamda34[i]*Q[i,3]*Q[i,4]
    lamda35[i]<-xlamda35[i]*Q[i,3]*Q[i,5]
    lamda45[i]<-xlamda45[i]*Q[i,4]*Q[i,5]
  }
  lamda<-cbind(xlamda0,lamda1,lamda2,lamda3,lamda4,lamda5,
               lamda12,lamda13,lamda14,lamda15,lamda23,lamda24,lamda25,lamda34,lamda35,lamda45
  )
  ### Probability of respondents in all attribute mastery patterns on all paired items
  P_KS <- matrix(NA,nrow = C,ncol = Npair)
  for (c in 1:C) {
    for (b in 1:Npair) {
      ta <- xlamda0[pair_in[b,1]]+xlamda1[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]+xlamda2[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]+xlamda3[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda4[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda5[pair_in[b,1]]*all.patterns[c,5]*Q[pair_in[b,1],5]
      +xlamda12[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,2]*Q[pair_in[b,1],2]+xlamda13[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda14[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda15[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,5]*Q[pair_in[b,1],5]+
        xlamda23[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda24[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda25[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda34[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,4]*Q[pair_in[b,1],4]+
        xlamda35[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda45[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]*all.patterns[c,5]*Q[pair_in[b,1],5]
      tb <- xlamda0[pair_in[b,2]]+xlamda1[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]+xlamda2[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]+xlamda3[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda4[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda5[pair_in[b,2]]*all.patterns[c,5]*Q[pair_in[b,2],5]
      +xlamda12[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,2]*Q[pair_in[b,2],2]+xlamda13[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda14[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda15[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,5]*Q[pair_in[b,2],5]+
        xlamda23[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda24[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda25[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda34[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,4]*Q[pair_in[b,2],4]+
        xlamda35[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda45[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]*all.patterns[c,5]*Q[pair_in[b,2],5]
      P_KS[c,b]<-exp(ta - tb)/(1+exp(ta - tb))
      
    }
  }
  
  ### Probability of all respondents on all paired items
  item_true<-cbind(all.patterns,P_KS)
  location<-NULL
  for (n in 1:Nperson) {
    temp<-rowSums(matrix(true_alpha[n,],nrow(all.patterns),ncol(all.patterns),byrow = T)==item_true[,1:K])
    location[n]<-which(temp==ncol(all.patterns))  
  }
  P <-item_true[location,-c(1:5)]
  Y <- matrix(NA,Nperson,Npair)
  ### generate response data
  for(n in 1:Nperson){
    for(b in 1:Npair){
      if(runif(1,min=0,max=1) < P[n,b]){
        Y[n,b] <- 1
      }else{
        Y[n,b] <- 0
      }
    }
  }
  
  # Save Simulation Data--------------------------------------------------------------
  data <- new.env()
  data$Nperson <- Nperson
  data$Nstatement <- Nstatement
  data$Npair <- Npair
  data$Nblock <- Nblock
  data$BS <- BS
  data$K <- K
  data$C <- C
  data$R <- R
  data$Q <- Q                        
  data$pair_in <- pair_in
  data$true_alpha <- true_alpha      
  data$Y <- as.matrix(Y)             
  data$all.patterns<-all.patterns
  data$P_KS<-item_true               
  data$lamda<-lamda                  
  data$item_mu <- item_mu            
  data$item_sigma<-item_sigma        
  data <- as.list(data)
  # save data
  dput(data, file = path)
  # load data - dget()
  return(data)
}
