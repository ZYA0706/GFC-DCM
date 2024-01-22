################################################################################
#                                                                              #
#     Generalized Forced-Choice Diagnostic Classification Model Based on       #
#                                                                              #
#                   Thurstone's Law of Comparative Judgement                   #
#                                                                              #
#                           Simulation Study                                   #
#                                                                              #
################################################################################

# Notes：
# Please create the corresponding folder first to ensure that the program can run normally.
# Simulation_study>>simdata
# Simulation_study>>est
# Simulation_study>>result
# Simulation_study>>functions

# Step:1 Preparation------------------------------------------------------------
working='GFC_DCM_code_osf/Simulation_study'
setwd(working)
rm(list = ls())
gc()
## packages
library(GDINA)
library(gtools)
library(nimble)
## load functions
source(file = 'functions/para_recovery.R')
source(file = 'functions/develop_block.R')
source(file = 'functions/sim_GFC_DCM_5.R')
source(file = 'functions/sim_GFC_DCM_7.R')
source(file = 'functions/Est_GFC_DCM_5.R')
source(file = 'functions/Est_GFC_DCM_7.R')


# Step 2: simulation design ----------------------------------------------------
conditions <- expand.grid(
  Npersons = c(500,1000),         # sample size : 300, 500, 1000
  Nstatement = c(30, 60),         # number of statements : 30, 60
  BS = c(2,3),                    # FC form : pair, triplet
  K = c(5, 7),                    # number of attributes : 5, 7
  LSM = c('unif','mvnorm'),       # latent structrue model : uniform, mvnorm
  stringsAsFactors = FALSE
)
conditions$sim_name = do.call("paste0", list("Nperson=", conditions$Npersons,
                                             "_Nstatement=",conditions$Nstatement,
                                             "_BS=", conditions$BS,
                                             "_K=",conditions$K,
                                             "_LSM=",conditions$LSM
))

# Step 3: simulate data and estimate parameters---------------------------------

Nrep<-30    #replication

time1 <- Sys.time()
for (s in 1:nrow(conditions)) {
  # result table
  PCCR<-matrix(NA,Nrep,1)
  colnames(PCCR)<-"Est_PCCR"
  if (conditions[s,]$K==5) {
    ACCR<-matrix(NA,Nrep,conditions[s,]$K)
    colnames(ACCR)<-c("Est_A1","Est_A2","Est_A3","Est_A4","Est_A5")
  }else {
    ACCR<-matrix(NA,Nrep,conditions[s,]$K)
    colnames(ACCR)<-c("Est_A1","Est_A2","Est_A3","Est_A4","Est_A5","Est_A6","Est_A7")
  }
  Esti_Intercept<-matrix(NA,Nrep,4)
  colnames(Esti_Intercept)<-c("Bias_Intercept","MAE_Intercept","MSE_Intercept","RMSE_Intercept")
  Esti_Main<-matrix(NA,Nrep,4)
  colnames(Esti_Main)<-c("Bias_Main","MAE_Main","MSE_Main","RMSE_Main")
  Esti_Interaction<-matrix(NA,Nrep,4)
  colnames(Esti_Interaction)<-c("Bias_Interaction","MAE_Interaction","MSE_Interaction","RMSE_Interaction")
  for (rep in 1:Nrep) {
    cond_K <- conditions[s,]$K
    if (cond_K == 5) {
      # simulation
      simdata <- sim_GFC_DCM_5(s = s,conditions = conditions,rep = rep)
      # estimation
      MCMCdata <- Est_GFC_DCM_5(simdata = simdata,
                                s = s,conditions = conditions,rep = rep)
    }else if(cond_K == 7){
      # simulation
      simdata <- sim_GFC_DCM_7(s = s,conditions = conditions,rep = rep)
      # estimation
      MCMCdata <- Est_GFC_DCM_7(simdata = simdata,
                                s = s,conditions = conditions,rep = rep)
    }
    PCCR[rep,] <- MCMCdata$PCCR
    ACCR[rep,] <- MCMCdata$ACCR
    Esti_Intercept[rep,] <- MCMCdata$Esti_Intercept
    Esti_Main[rep,] <- MCMCdata$Esti_Main
    Esti_Interaction[rep,] <- MCMCdata$Esti_Interaction
    gc()
    print(PCCR[rep,])
    print(ACCR[rep,])
    print(Esti_Intercept[rep,])
    print(Esti_Main[rep,])
    print(Esti_Interaction[rep,])
  }
  # save results
  write.csv(PCCR, file = paste("result/",'PCCR_',conditions$sim_name[s],'.csv',sep=' '))
  write.csv(ACCR, file = paste("result/",'ACCR_',conditions$sim_name[s],'.csv',sep=' '))
  write.csv(Esti_Intercept, file = paste("result/",'Esti_Intercept_',conditions$sim_name[s],'.csv',sep=' '))
  write.csv(Esti_Main, file = paste("result/",'Esti_Main',conditions$sim_name[s],'.csv',sep=' '))
  write.csv(Esti_Interaction, file = paste("result/",'Esti_Interaction',conditions$sim_name[s],'.csv',sep=' '))
}

time2 <- Sys.time()
used_time <- time2 - time1

## end

