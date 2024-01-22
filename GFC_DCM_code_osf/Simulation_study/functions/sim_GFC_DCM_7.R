# This function was used to generate the data when K=7 in simulation.
sim_FC_TDCM_7 <- function(s,conditions,rep){
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
  R <- diag(29)
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
  xlamda6_mu <- 4
  xlamda7_mu <- 4
  xlamda12_mu <- 0
  xlamda13_mu <- 0
  xlamda14_mu <- 0
  xlamda15_mu <- 0
  xlamda16_mu <- 0
  xlamda17_mu <- 0
  xlamda23_mu <- 0
  xlamda24_mu <- 0
  xlamda25_mu <- 0
  xlamda26_mu <- 0
  xlamda27_mu <- 0
  xlamda34_mu <- 0
  xlamda35_mu <- 0
  xlamda36_mu <- 0
  xlamda37_mu <- 0
  xlamda45_mu <- 0
  xlamda46_mu <- 0
  xlamda47_mu <- 0
  xlamda56_mu <- 0
  xlamda57_mu <- 0
  xlamda67_mu <- 0
  
  
  item_mu <- c(xlamda0_mu,xlamda1_mu,xlamda2_mu,xlamda3_mu,xlamda4_mu,xlamda5_mu,xlamda6_mu,xlamda7_mu,
               xlamda12_mu,xlamda13_mu,xlamda14_mu,xlamda15_mu,xlamda16_mu,xlamda17_mu,
               xlamda23_mu,xlamda24_mu,xlamda25_mu,xlamda26_mu,xlamda27_mu,
               xlamda34_mu,xlamda35_mu,xlamda36_mu,xlamda37_mu,
               xlamda45_mu,xlamda46_mu,xlamda47_mu,
               xlamda56_mu,xlamda57_mu,
               xlamda67_mu)
  
  ### covariance matrix
  
  item_sigma <- diag(29)
  diag(item_sigma) <- 0.5
  item_mvrnorm <- MASS::mvrnorm(Nstatement, mu = item_mu, Sigma = item_sigma, 29)
  
  xlamda0 <- as.matrix(item_mvrnorm[,1])
  xlamda1 <- as.matrix(item_mvrnorm[,2])
  xlamda2 <- as.matrix(item_mvrnorm[,3])
  xlamda3 <- as.matrix(item_mvrnorm[,4])
  xlamda4 <- as.matrix(item_mvrnorm[,5])
  xlamda5 <- as.matrix(item_mvrnorm[,6])
  xlamda6 <- as.matrix(item_mvrnorm[,7])
  xlamda7 <- as.matrix(item_mvrnorm[,8])
  xlamda12 <- as.matrix(item_mvrnorm[,9])
  xlamda13 <- as.matrix(item_mvrnorm[,10])
  xlamda14 <- as.matrix(item_mvrnorm[,11])
  xlamda15 <- as.matrix(item_mvrnorm[,12])
  xlamda16 <- as.matrix(item_mvrnorm[,13])
  xlamda17 <- as.matrix(item_mvrnorm[,14])
  xlamda23 <- as.matrix(item_mvrnorm[,15])
  xlamda24 <- as.matrix(item_mvrnorm[,16])
  xlamda25 <- as.matrix(item_mvrnorm[,17])
  xlamda26 <- as.matrix(item_mvrnorm[,18])
  xlamda27 <- as.matrix(item_mvrnorm[,19])
  xlamda34 <- as.matrix(item_mvrnorm[,20])
  xlamda35 <- as.matrix(item_mvrnorm[,21])
  xlamda36 <- as.matrix(item_mvrnorm[,22])
  xlamda37 <- as.matrix(item_mvrnorm[,23])
  xlamda45 <- as.matrix(item_mvrnorm[,24])
  xlamda46 <- as.matrix(item_mvrnorm[,25])
  xlamda47 <- as.matrix(item_mvrnorm[,26])
  xlamda56 <- as.matrix(item_mvrnorm[,27])
  xlamda57 <- as.matrix(item_mvrnorm[,28])
  xlamda67 <- as.matrix(item_mvrnorm[,29])
  
  lamda1<-matrix(NA,Nstatement,1)
  lamda2<-matrix(NA,Nstatement,1)
  lamda3<-matrix(NA,Nstatement,1)
  lamda4<-matrix(NA,Nstatement,1)
  lamda5<-matrix(NA,Nstatement,1)
  lamda6<-matrix(NA,Nstatement,1)
  lamda7<-matrix(NA,Nstatement,1)
  lamda12<-matrix(NA,Nstatement,1)
  lamda13<-matrix(NA,Nstatement,1)
  lamda14<-matrix(NA,Nstatement,1)
  lamda15<-matrix(NA,Nstatement,1)
  lamda16<-matrix(NA,Nstatement,1)
  lamda17<-matrix(NA,Nstatement,1)
  lamda23<-matrix(NA,Nstatement,1)
  lamda24<-matrix(NA,Nstatement,1)
  lamda25<-matrix(NA,Nstatement,1)
  lamda26<-matrix(NA,Nstatement,1)
  lamda27<-matrix(NA,Nstatement,1)
  lamda34<-matrix(NA,Nstatement,1)
  lamda35<-matrix(NA,Nstatement,1)
  lamda36<-matrix(NA,Nstatement,1)
  lamda37<-matrix(NA,Nstatement,1)
  lamda45<-matrix(NA,Nstatement,1)
  lamda46<-matrix(NA,Nstatement,1)
  lamda47<-matrix(NA,Nstatement,1)
  lamda56<-matrix(NA,Nstatement,1)
  lamda57<-matrix(NA,Nstatement,1)
  lamda67<-matrix(NA,Nstatement,1)
  
  for(i in 1:Nstatement){
    lamda1[i]<-xlamda1[i]*Q[i,1]
    lamda2[i]<-xlamda2[i]*Q[i,2]
    lamda3[i]<-xlamda3[i]*Q[i,3]
    lamda4[i]<-xlamda4[i]*Q[i,4]
    lamda5[i]<-xlamda5[i]*Q[i,5]
    lamda6[i]<-xlamda6[i]*Q[i,6]
    lamda7[i]<-xlamda7[i]*Q[i,7]
    lamda12[i]<-xlamda12[i]*Q[i,1]*Q[i,2]
    lamda13[i]<-xlamda13[i]*Q[i,1]*Q[i,3]
    lamda14[i]<-xlamda14[i]*Q[i,1]*Q[i,4]
    lamda15[i]<-xlamda15[i]*Q[i,1]*Q[i,5]
    lamda16[i]<-xlamda16[i]*Q[i,1]*Q[i,6]
    lamda17[i]<-xlamda17[i]*Q[i,1]*Q[i,7]
    lamda23[i]<-xlamda23[i]*Q[i,2]*Q[i,3]
    lamda24[i]<-xlamda24[i]*Q[i,2]*Q[i,4]
    lamda25[i]<-xlamda25[i]*Q[i,2]*Q[i,5]
    lamda26[i]<-xlamda26[i]*Q[i,2]*Q[i,6]
    lamda27[i]<-xlamda27[i]*Q[i,2]*Q[i,7]
    lamda34[i]<-xlamda34[i]*Q[i,3]*Q[i,4]
    lamda35[i]<-xlamda35[i]*Q[i,3]*Q[i,5]
    lamda36[i]<-xlamda36[i]*Q[i,3]*Q[i,6]
    lamda37[i]<-xlamda37[i]*Q[i,3]*Q[i,7]
    lamda45[i]<-xlamda45[i]*Q[i,4]*Q[i,5]
    lamda46[i]<-xlamda46[i]*Q[i,4]*Q[i,6]
    lamda47[i]<-xlamda47[i]*Q[i,4]*Q[i,7]
    lamda56[i]<-xlamda56[i]*Q[i,5]*Q[i,6]
    lamda57[i]<-xlamda57[i]*Q[i,5]*Q[i,7]
    lamda67[i]<-xlamda67[i]*Q[i,6]*Q[i,7]
  }
  lamda<-cbind(xlamda0,lamda1,lamda2,lamda3,lamda4,lamda5,lamda6,lamda7,
               lamda12,lamda13,lamda14,lamda15,lamda16,lamda17,
               lamda23,lamda24,lamda25,lamda26,lamda27,
               lamda34,lamda35,lamda36,lamda37,
               lamda45,lamda46,lamda47,
               lamda56,lamda57,
               lamda67
  )
  ### Probability of respondents in all attribute mastery patterns on all paired items
  P_KS <- matrix(NA,nrow = C,ncol = Npair)
  for (c in 1:C) {
    for (b in 1:Npair) {
      ta <- xlamda0[pair_in[b,1]]+xlamda1[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]+xlamda2[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]+xlamda3[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda4[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda5[pair_in[b,1]]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda6[pair_in[b,1]]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda7[pair_in[b,1]]*all.patterns[c,7]*Q[pair_in[b,1],7]
      +xlamda12[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,2]*Q[pair_in[b,1],2]+xlamda13[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda14[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda15[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda16[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda17[pair_in[b,1]]*all.patterns[c,1]*Q[pair_in[b,1],1]*all.patterns[c,7]*Q[pair_in[b,1],7]+
        xlamda23[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,3]*Q[pair_in[b,1],3]+xlamda24[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,4]*Q[pair_in[b,1],4]+xlamda25[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda26[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda27[pair_in[b,1]]*all.patterns[c,2]*Q[pair_in[b,1],2]*all.patterns[c,7]*Q[pair_in[b,1],7]+xlamda34[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,4]*Q[pair_in[b,1],4]+
        xlamda35[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda36[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda37[pair_in[b,1]]*all.patterns[c,3]*Q[pair_in[b,1],3]*all.patterns[c,7]*Q[pair_in[b,1],7]+xlamda45[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]*all.patterns[c,5]*Q[pair_in[b,1],5]+xlamda46[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda47[pair_in[b,1]]*all.patterns[c,4]*Q[pair_in[b,1],4]*all.patterns[c,7]*Q[pair_in[b,1],7]+xlamda56[pair_in[b,1]]*all.patterns[c,5]*Q[pair_in[b,1],5]*all.patterns[c,6]*Q[pair_in[b,1],6]+xlamda57[pair_in[b,1]]*all.patterns[c,5]*Q[pair_in[b,1],5]*all.patterns[c,7]*Q[pair_in[b,1],7]+xlamda67[pair_in[b,1]]*all.patterns[c,6]*Q[pair_in[b,1],6]*all.patterns[c,7]*Q[pair_in[b,1],7]
      tb <- xlamda0[pair_in[b,2]]+xlamda1[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]+xlamda2[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]+xlamda3[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda4[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda5[pair_in[b,2]]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda6[pair_in[b,2]]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda7[pair_in[b,2]]*all.patterns[c,7]*Q[pair_in[b,2],7]
      +xlamda12[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,2]*Q[pair_in[b,2],2]+xlamda13[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda14[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda15[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda16[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda17[pair_in[b,2]]*all.patterns[c,1]*Q[pair_in[b,2],1]*all.patterns[c,7]*Q[pair_in[b,2],7]+
        xlamda23[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,3]*Q[pair_in[b,2],3]+xlamda24[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,4]*Q[pair_in[b,2],4]+xlamda25[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda26[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda27[pair_in[b,2]]*all.patterns[c,2]*Q[pair_in[b,2],2]*all.patterns[c,7]*Q[pair_in[b,2],7]+xlamda34[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,4]*Q[pair_in[b,2],4]+
        xlamda35[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda36[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda37[pair_in[b,2]]*all.patterns[c,3]*Q[pair_in[b,2],3]*all.patterns[c,7]*Q[pair_in[b,2],7]+xlamda45[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]*all.patterns[c,5]*Q[pair_in[b,2],5]+xlamda46[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda47[pair_in[b,2]]*all.patterns[c,4]*Q[pair_in[b,2],4]*all.patterns[c,7]*Q[pair_in[b,2],7]+xlamda56[pair_in[b,2]]*all.patterns[c,5]*Q[pair_in[b,2],5]*all.patterns[c,6]*Q[pair_in[b,2],6]+xlamda57[pair_in[b,2]]*all.patterns[c,5]*Q[pair_in[b,2],5]*all.patterns[c,7]*Q[pair_in[b,2],7]+xlamda67[pair_in[b,2]]*all.patterns[c,6]*Q[pair_in[b,2],6]*all.patterns[c,7]*Q[pair_in[b,2],7]
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
  P <-item_true[location,-c(1:7)]
  ### generate response data
  Y <- matrix(NA,Nperson,Npair)
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
