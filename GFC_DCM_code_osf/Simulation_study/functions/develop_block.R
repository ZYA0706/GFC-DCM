develop_block2345<-function(item_number,Block_size = 2){
  library(gtools)
  allblock<-t(combn(item_number,Block_size))
  final_block<-matrix(NA,nrow=length(item_number)/Block_size,ncol=Block_size)
  if (Block_size == 2) {
    for(k in 1:nrow(final_block)){
      blockchu<-allblock[sample(nrow(allblock),1,replace=F),]#??????block
      final_block[k,]<-as.vector(blockchu)
      truehang<-which(allblock[,1]==blockchu[1]|allblock[,1]==blockchu[2]|
                        allblock[,2]==blockchu[1]|allblock[,2]==blockchu[2])
      remain_block<-as.matrix(allblock[-c(truehang),])
      if(length(remain_block)==2){
        final_block[nrow(final_block),]<-remain_block
        break
      }else{
        allblock<-remain_block
      }
    }
  }else if(Block_size == 3){
    for(k in 1:nrow(final_block)){
      blockchu<-allblock[sample(nrow(allblock),1,replace=F),]#??????block
      final_block[k,]<-as.vector(blockchu)
      truehang<-which(allblock[,1]==blockchu[1]|allblock[,1]==blockchu[2]|allblock[,1]==blockchu[3]|
                        allblock[,2]==blockchu[1]|allblock[,2]==blockchu[2]|allblock[,2]==blockchu[3]|
                        allblock[,3]==blockchu[1]|allblock[,3]==blockchu[2]|allblock[,3]==blockchu[3]
      )
      remain_block<-as.matrix(allblock[-c(truehang),])
      if(length(remain_block)==3){
        final_block[nrow(final_block),]<-remain_block
        break
      }else{
        allblock<-remain_block
      }
    }
  }else if(Block_size == 4){
    for(k in 1:nrow(final_block)){
      blockchu<-allblock[sample(nrow(allblock),1,replace=F),]#??????block
      final_block[k,]<-as.vector(blockchu)
      truehang<-which(allblock[,1]==blockchu[1]|allblock[,1]==blockchu[2]|allblock[,1]==blockchu[3]|allblock[,1]==blockchu[4]|
                        allblock[,2]==blockchu[1]|allblock[,2]==blockchu[2]|allblock[,2]==blockchu[3]|allblock[,2]==blockchu[4]|
                        allblock[,3]==blockchu[1]|allblock[,3]==blockchu[2]|allblock[,3]==blockchu[3]|allblock[,3]==blockchu[4]|
                        allblock[,4]==blockchu[1]|allblock[,4]==blockchu[2]|allblock[,4]==blockchu[3]|allblock[,4]==blockchu[4])
      remain_block<-as.matrix(allblock[-c(truehang),])
      if(length(remain_block)==4){
        final_block[nrow(final_block),]<-remain_block
        break
      }else{
        allblock<-remain_block
      }
    }
  }else if(Block_size == 5){
    for(k in 1:nrow(final_block)){
      blockchu<-allblock[sample(nrow(allblock),1,replace=F),]#??????block
      final_block[k,]<-as.vector(blockchu)
      truehang<-which(allblock[,1]==blockchu[1]|allblock[,1]==blockchu[2]|allblock[,1]==blockchu[3]|allblock[,1]==blockchu[4]|allblock[,1]==blockchu[5]|
                        allblock[,2]==blockchu[1]|allblock[,2]==blockchu[2]|allblock[,2]==blockchu[3]|allblock[,2]==blockchu[4]|allblock[,2]==blockchu[5]|
                        allblock[,3]==blockchu[1]|allblock[,3]==blockchu[2]|allblock[,3]==blockchu[3]|allblock[,3]==blockchu[4]|allblock[,3]==blockchu[5]|
                        allblock[,4]==blockchu[1]|allblock[,4]==blockchu[2]|allblock[,4]==blockchu[3]|allblock[,4]==blockchu[4]|allblock[,4]==blockchu[5]|
                        allblock[,5]==blockchu[1]|allblock[,5]==blockchu[2]|allblock[,5]==blockchu[3]|allblock[,5]==blockchu[4]|allblock[,5]==blockchu[5]
      )
      remain_block<-as.matrix(allblock[-c(truehang),])
      if(length(remain_block)==5){
        final_block[nrow(final_block),]<-remain_block
        break
      }else{
        allblock<-remain_block
      }
    }
  }
  
  
  return(final_block)
}
