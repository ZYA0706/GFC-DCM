################################################################################
#                                                                              #
#     Generalized Forced-Choice Diagnostic Classification Model Based on       #
#                                                                              #
#                   Thurstone's Law of Comparative Judgement                   #
#                                                                              #
#                       Empirical Study Result Summary                         #
#                                                                              #
################################################################################

# Notes：
# Please create the corresponding folder first to ensure that the program can run normally.
# Empirical_study>>realdata
# Empirical_study>>result

# Step:1 Preparation------------------------------------------------------------
working='GFC_DCM_code_osf/Empirical study'
setwd(working)
rm(list = ls())
gc()

## load packages
library(ggplot2)
library(GDINA)

## load result lists
GFC_DCM <- readRDS(file = 'result/RD_GFC_DCM_list.rds')
FC_DCM <- readRDS(file = 'result/RD_FC_DCM_list.rds')


## nessessary infomation about realdata
Nperson <- GFC_DCM$Nperson;
Nstatement <- GFC_DCM$Nstatement;
Nblock <- GFC_DCM$Nblock
Npair <- GFC_DCM$Npair;
K <- GFC_DCM$K
C <- 2^K
Q <- GFC_DCM$Q;
all.patterns <- GFC_DCM$data$all.patterns
block <- matrix(1:Nstatement,nrow = Nblock,byrow = T)
pair_in <- NULL
for (i in 1:Nblock) {
  pair_in <- append(x = pair_in,as.vector(combn(x = block[i,],2)))
}
pair_in <- matrix(pair_in,ncol = 2,byrow = T)
Q1 = Q[pair_in[,1],]
Q2 = Q[pair_in[,2],]
Y <- GFC_DCM$data$Y   # paired items response data

# Step:2 Statistical Analysis---------------------------------------------------
## WAIC-------------------------------------------------------------------------

GFC_DCM_WAIC <- GFC_DCM$WAIC;GFC_DCM_WAIC$WAIC
FC_DCM_WAIC <- FC_DCM$WAIC;FC_DCM_WAIC$WAIC

## PPP-value--------------------------------------------------------------------

### test-level PPP-value
GFC_DCM$ppp_test;FC_DCM$ppp_test

### block-level PPP-value
GFC_DCM_PPP <- round(GFC_DCM$ppp_block,3);GFC_DCM_PPP
FC_DCM_PPP <- round(FC_DCM$ppp_block,3);FC_DCM_PPP

### ploting block-level PPP-value for both models
df <- cbind(NO_BLOCK = rep(1:11,2),
            Model = rep(c('GFC-DCM',"FC-DCN"),each = 11),
            PPP = c(GFC_DCM_PPP,FC_DCM_PPP))
df <- as.data.frame(df)
df$Model <- factor(df$Model, levels = c('GFC-DCM',"FC-DCN"))
df$NO_BLOCK<- as.numeric(df$NO_BLOCK)
df$PPP <- as.numeric(df$PPP)
dev.new()
p <- ggplot(df, aes(x = NO_BLOCK, y = PPP, color = Model)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black",size = 0.6) +
  # theme_void() +  # 设置图表背景为空白
  # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), sec.axis = sec_axis(~., breaks = seq(0, 1, 0.1))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 11, 1)) + 
  labs(x = "No.Block", y = "PPP") +
  scale_color_manual(values = c("green","red"),labels = c('GFC-DCM',"FC-DCN")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 14,face = 'bold'),
        axis.title = element_text(size = 16,face = 'bold'),
        legend.text = element_text(size = 14),  # 调整图例文本的大小
        legend.title = element_text(size = 16))  # 调整图例标题的大小  # 将图例放在底部 调整刻度数字的大小

p
ggsave("result/PPP-value_for_two_models.png", plot = p, device = "png",bg = 'white',width = 12, height = 6, units = 'in')

## Attribute reliability--------------------------------------------------------

alpha_GFC_DCM <- GFC_DCM$est_alpha
alpha_FC_DCM <- FC_DCM$est_alpha

### GFC_DCM---------------------------------------------------------------------

### mixed proportion for attribute mastery pattern

classify_matrix <- function(alpha,all.patterns=all.patterns){
  N <- nrow(alpha)
  K <- ncol(all.patterns)
  C <- 2^K
  class <- NULL
  for (n in 1:N) {
    index <- matrix(NA,C,K)
    for (c in 1:C) {
      index[c,] <- alpha[n,] == all.patterns[c,]
    }
    class[n] <- which(apply(index,1,all))
  }
  return(class)
}
est_class <- classify_matrix(alpha = alpha_GFC_DCM,all.patterns = all.patterns)
mixed_proportion <- table(est_class)/Nperson


### success probability of all attribute mastery patterns
P_c <- matrix(NA,C,Npair)
w1 <- w2 <- array(NA,dim = c(C,Npair,4))
alpha <- GFC_DCM$data$all.patterns
lamda <- GFC_DCM$est_lamda
s_eta1 <- s_eta2 <- s_eta3 <- k_eta1 <- k_eta2 <- k_eta3 <- matrix(NA,C,Npair)

for(c in 1:C){
  for(b in 1:Npair){
    for (k in 1:K){
      w1[c, b, k] <-alpha[c, k]^Q1[b, k]
      w2[c, b, k] <- alpha[c, k]^Q2[b, k]}
    # first statemect in paired item
    s_eta1[c, b] <- lamda[,2][pair_in[b,1]] * w1[c, b, 1] + lamda[,3][pair_in[b,1]] * w1[c, b, 2] + lamda[,4][pair_in[b,1]] * w1[c, b, 3] + lamda[,5][pair_in[b,1]] * w1[c, b, 4]
    s_eta2[c, b] <- lamda[,6][pair_in[b,1]] * w1[c, b, 1] * w1[c, b, 2] + lamda[,7][pair_in[b,1]] * w1[c, b, 1] * w1[c, b, 3] + lamda[,8][pair_in[b,1]] * w1[c, b, 1] * w1[c, b, 4] + lamda[,9][pair_in[b,1]] * w1[c, b, 2] * w1[c, b, 4] + lamda[,10][pair_in[b,1]] * w1[c, b, 3] * w1[c, b, 4]
    s_eta3[c, b] <- lamda[,11][pair_in[b,1]] * w1[c, b, 2] * w1[c, b, 3] * w1[c, b, 4]
    
    # secocd statemect ic paired item
    k_eta1[c, b] <- lamda[,2][pair_in[b,2]] * w2[c, b, 1] + lamda[,3][pair_in[b,2]] * w2[c, b, 2] + lamda[,4][pair_in[b,2]] * w2[c, b, 3] + lamda[,5][pair_in[b,2]] * w2[c, b, 4]
    k_eta2[c, b] <- lamda[,6][pair_in[b,2]] * w2[c, b, 1] * w2[c, b, 2] + lamda[,7][pair_in[b,2]] * w2[c, b, 1] * w2[c, b, 3]  + lamda[,8][pair_in[b,2]] * w2[c, b, 1] * w2[c, b, 4] + lamda[,9][pair_in[b,2]] * w2[c, b, 2] * w2[c, b, 4] + lamda[,10][pair_in[b,2]] * w2[c, b, 3] * w2[c, b, 4]
    k_eta3[c, b] <- lamda[,11][pair_in[b,2]] * w2[c, b, 2] * w2[c, b, 3] * w2[c, b, 4]
    
    
    P_c[c,b] <- nimble::ilogit((lamda[,1][pair_in[b,1]] + s_eta1[c, b] + s_eta2[c, b] + s_eta3[c, b]) - (lamda[,1][pair_in[b,2]] + k_eta1[c, b] + k_eta2[c, b] + k_eta3[c, b]))
  }
}

P <-  liklihood <- matrix(NA,Nperson,Npair)
for (n in 1:Nperson) {
  liklihood[n,] <- P[n,] <- P_c[est_class[n],]
  for (b in 1:Npair) {
    if (Y[n,b]==0) {
      liklihood[n,b] <- (1 - P[n,b])
    }
  }
}


## respondent attribute mastery pattern probability estimate

aic <- matrix(NA,nrow = Nperson,ncol = C)

for (n in 1:Nperson) {
  dom <- NULL
  for (c in 1:C) {
    dom[c] <- mixed_proportion[c]*prod((P_c[c,]^Y[n,] )* ((1-P_c[c,])^(1-Y[n,])))
  }
  aic[n,] <- dom/sum(dom)
}


## Respondent-level marginal probability of attribute mastery


pik_TDCM <- matrix(NA,nrow = Nperson, K)

for (i in 1:Nperson) {
  pik_TDCM[i,] <- t(GFC_DCM$data$all.patterns) %*% as.matrix(aic[i,])
}



#### computer attribute reliability-----------------------------------------------

r_TDCM <- a <- b <- c <- d <-  NULL
for (k in 1:K) {
  pp <- pik_TDCM[,k] 
  a[k] <- mean(pp*pp) * Nperson
  d[k] <- mean((1-pp)*(1-pp)) * Nperson
  b[k] <- c[k] <- mean(pp*(1-pp)) * Nperson
  r_TDCM[k] <- cos(pi/(1 + sqrt((a[k]*d[k])/(b[k]*c[k]))))
}


### FC-DCM----------------------------------------------------------------------

#### Respondent-level marginal probability of attribute mastery-------------------

pik_DCM <- matrix(NA,nrow = Nperson, K)
theta <- FC_DCM$est_theta
delta0 <- FC_DCM$est_delta0
delta1 <- FC_DCM$est_delta1
for (n in 1:Nperson) {
  for (k in 1:K) {
    pik_DCM[n,k] <- nimble::ilogit(1.7*delta1[k]*(theta[n]-delta0[k]))
  }
}


#### computer attribute reliability-----------------------------------------------

r_DCM <- a <- b <- c <- d <-  NULL
for (k in 1:K) {
  pp <- pik_DCM[,k] 
  a[k] <- (mean(pp*pp) * Nperson)
  d[k] <- (mean((1-pp)*(1-pp)) * Nperson)
  b[k] <- c[k] <- (mean(pp*(1-pp)) * Nperson)
  r_DCM[k] <- cos(pi/(1 + sqrt((a[k]*d[k])/(b[k]*c[k]))))
}

### ploting attribute-level reliability for both models
df <- cbind(NO_Attribute = rep(1:4,2),
            Model = rep(c('GFC-DCM',"FC-DCN"),each = 4),
            Rel = c(r_TDCM,r_DCM))
df <- as.data.frame(df)
df$Model <- factor(df$Model, levels = c('GFC-DCM',"FC-DCN"))
df$NO_Attribute<- as.numeric(df$NO_Attribute)
df$Rel <- as.numeric(df$Rel)
dev.new()
p <- ggplot(df, aes(x = NO_Attribute, y = Rel, color = Model)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black",size = 0.6) +
  # theme_void() +  # 设置图表背景为空白
  # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), sec.axis = sec_axis(~., breaks = seq(0, 1, 0.1))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 4, 1)) + 
  labs(x = "No.Attribute", y = "Reliability") +
  scale_color_manual(values = c("green","red"),labels = c('GFC-DCM',"FC-DCN")) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text = element_text(size = 14,face = 'bold'),
        axis.title = element_text(size = 16,face = 'bold'),
        legend.text = element_text(size = 14),  # 调整图例文本的大小
        legend.title = element_text(size = 16))  # 调整图例标题的大小  # 将图例放在底部 调整刻度数字的大小

p
ggsave("Rel-value_for_two_models.png", plot = p, device = "png",bg = 'white',width = 12, height = 6, units = 'in')



## Esitimated Parameter Analysis------------------------------------------------

### Attribute prevalence result both models-------------------------------------
apply(alpha_FC_DCM, 2, sum)/Nperson
apply(alpha_GFC_DCM, 2, sum)/Nperson
colMeans(pik_DCM)
colMeans(pik_TDCM)
### agreement ratios of estimated attribute score results both models-----------
sum(alpha_FC_DCM[,1] == alpha_GFC_DCM[,1])/Nperson
sum(alpha_FC_DCM[,2] == alpha_GFC_DCM[,2])/Nperson
sum(alpha_FC_DCM[,3] == alpha_GFC_DCM[,3])/Nperson
sum(alpha_FC_DCM[,4] == alpha_GFC_DCM[,4])/Nperson

### correlations between the marginal mastery probability of attribute and estimated attribute score----------

#### GFC_DCM
cor.test(pik_TDCM[,1],alpha_GFC_DCM[,1])
cor.test(pik_TDCM[,2],alpha_GFC_DCM[,2])
cor.test(pik_TDCM[,3],alpha_GFC_DCM[,3])
cor.test(pik_TDCM[,4],alpha_GFC_DCM[,4])

#### FC-DCM
cor.test(pik_DCM[,1],alpha_FC_DCM[,1])
cor.test(pik_DCM[,2],alpha_FC_DCM[,2])
cor.test(pik_DCM[,3],alpha_FC_DCM[,3])
cor.test(pik_DCM[,4],alpha_FC_DCM[,4])

### statements parameters-------------------------------------------------------

#### GFC_DCM
intercept_effect <- GFC_DCM$est_lamda[,1]
range(intercept_effect);mean(intercept_effect)
write.csv(matrix(intercept_effect,Nblock,byrow = T,
                 dimnames = list( paste0('block',1:Nblock))),
          file = 'result/intercept_effect.csv')
write.csv(matrix(GFC_DCM$psd_lamda[,1],Nblock,byrow = T,
                 dimnames = list( paste0('block',1:Nblock))),
          file = 'result/intercept_effect_SE.csv')
mean(GFC_DCM$psd_lamda[,1]) # mean psd of intercept effect


main_effect <- GFC_DCM$est_lamda[,2:5]
main_effect[main_effect==0] <- NA
rownames(main_effect) <- rep(paste0('block',1:Nblock),each=4)
apply(main_effect, 2, function(x){range(x,na.rm = T)});apply(main_effect,2,function(x){mean.default(x,na.rm = T)})
write.csv(main_effect,file = 'result/main_effect.csv')
write.csv(GFC_DCM$psd_lamda[,2:5],file = 'result/main_effect_SE.csv')
sum(GFC_DCM$psd_lamda[,2:5])/61 # mean psd of main effect


interaction_effect <- GFC_DCM$est_lamda[,6:11]
interaction_effect[interaction_effect==0] <- NA
rownames(interaction_effect) <- rep(paste0('block',1:Nblock),each=4)
apply(interaction_effect, 2, function(x){range(x,na.rm = T)});apply(interaction_effect,2,function(x){mean.default(x,na.rm = T)})
range(interaction_effect,na.rm = T);mean.default(interaction_effect,na.rm = T)
write.csv(interaction_effect,file = 'result/interaction_effect.csv')
write.csv(GFC_DCM$psd_lamda[,6:11],
          file = 'result/interaction_effect_SE.csv')
sum(GFC_DCM$psd_lamda[,6:11])/18 # mean psd of interaction effect


#### FC-DCM

d0 <- FC_DCM$est_d0
write.csv(matrix(d0,Nblock,byrow = T),'result/d0.csv')
write.csv(matrix(FC_DCM$psd_d0,Nblock,byrow = T),'result/d0_SE.csv')
range(d0);mean(d0)
length(d0[d0>=0.2])/Npair # 0.8787879
length(d0[d0>=0.3])/Npair # 0.4848485
length(d0[d0>=0.4])/Npair # 0.2727273
mean(FC_DCM$psd_d0)


sd <- FC_DCM$est_sd
write.csv(matrix(sd,Nblock,byrow = T),'result/sd.csv')
write.csv(matrix(FC_DCM$psd_sd,Nblock,byrow = T),'result/sd_SE.csv')
range(sd);mean(sd)
mean(FC_DCM$psd_sd)


delta0 <- FC_DCM$est_delta0 
delta0
mean(FC_DCM$psd_delta0)


delta1 <- FC_DCM$est_delta1
delta1
mean(FC_DCM$psd_delta1)


theta <- FC_DCM$est_theta
range(theta);mean(theta)
mean(FC_DCM$psd_theta)
empirical_reliability_theta <- var(theta)/(var(theta)+mean(FC_DCM$psd_theta^2))


