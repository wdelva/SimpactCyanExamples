IDs.Seq.Random.skew <- function(simpact.trans.net = simpact.trans.net,
                                limitTransmEvents = 7,
                                timewindow = c(10,40),
                                seq.cov = 50,
                                age.limit=100,
                                age.group = c(25, 40),
                                propor=0.7){
  
  seeds.id <- length(simpact.trans.net)
  
  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)
  
  for(i in 1:seeds.id){
    
    transm.age.i <- as.data.frame(simpact.trans.net[[i]])
    
    age.i <- transm.age.i$SampTime - transm.age.i$TOBRec
    
    transm.age.i <- cbind(transm.age.i, age.i)
    
    new.transm.tab[[i]] <- transm.age.i
    
  }
  
  # ID numbers of Selected networks with at least limitTransmEvents + 1 indiviuals
  
  IDs.transm <- vector()
  
  TransmEventsCountVector <- vector()
  
  for(k in 1:seeds.id){
    trans.net.i.check <- as.data.frame(new.transm.tab[[k]])
    
    if(nrow(trans.net.i.check)>=limitTransmEvents){
      
      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
      
      IDs.transm <- c(IDs.transm, k)
    }
  }
  
  if(length(IDs.transm)>=1){
    
    ## Binding together all selected transmission transmission networks ##
    
    for (q in 1:length(IDs.transm)){
      
      if(q==1){
        p <- IDs.transm[q]
        trans.sum <- new.transm.tab[[p]]
        rename.id <- paste0(p,".",trans.sum$id,".C")
        trans.sum$id <- rename.id
        trans.sum.rename.id <- trans.sum
      }
      else{
        
        p <- IDs.transm[q]
        
        read.trans.sum <- new.transm.tab[[p]]
        rename.id.read <- paste0(p,".",read.trans.sum$id,".C")
        read.trans.sum$id <- rename.id.read
        trans.sum.rename.id <- rbind(trans.sum.rename.id, read.trans.sum)
      }
      
    }
    
    trans.sum.age.limit <- dplyr::filter(trans.sum.rename.id, trans.sum.rename.id$age.i<=age.limit)
    
    trans.sum.men <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="0" & trans.sum.age.limit$SampTime >= timewindow[1] & trans.sum.age.limit$SampTime <= timewindow[2])
    
    trans.sum.women <- dplyr::filter(trans.sum.age.limit, trans.sum.age.limit$GenderRec=="1" & trans.sum.age.limit$SampTime >= timewindow[1] & trans.sum.age.limit$SampTime <= timewindow[2])
    
    perc.100 <- nrow(trans.sum.men) + nrow(trans.sum.women) # total number of individuals with age limit
    
    trans.sum.men.women <- rbind(trans.sum.men, trans.sum.women)
    
    perc.seq.coverage <- round(perc.100*seq.cov/100) # total number of individuals at seq.cov sequence coverage
    
    ##
    
    perc.seq.coverage.skew <- round(perc.seq.coverage * propor) # number of individuals in the overweighted age group
    
    rem.skew <- perc.seq.coverage - perc.seq.coverage.skew # remaining individuals in other age groups
    
    
    # Men and women in the overweighted age group
    
    trans.sum.men.skew <- dplyr::filter(trans.sum.men, trans.sum.men$age.i>=25 & trans.sum.men$age.i < 40)
    
    trans.sum.women.skew <- dplyr::filter(trans.sum.women, trans.sum.women$age.i>=25 & trans.sum.women$age.i < 40)
    
    # 100 % in the overweighted age group
    perc.100.skew <- nrow(trans.sum.men.skew) + nrow(trans.sum.women.skew) # total number of individuals with age limit
    
    #
    trans.sum.men.women.skew <- rbind(trans.sum.men.skew, trans.sum.women.skew)
    
    samp.skew <- sample(trans.sum.men.women.skew$id, perc.seq.coverage.skew) # sample in the overweighted age group
    
    
    # Data table of remaing age groups
    
    id.diff <- setdiff(trans.sum.men.women$id, trans.sum.men.women.skew$id)
    
    trans.sum.men.women.rem <- dplyr::filter(trans.sum.men.women, trans.sum.men.women$id %in%id.diff)
    
    samp.rem <- sample(trans.sum.men.women.rem$id, rem.skew) # sample in the remaining age groups
    
    samp.all <- c(samp.skew, samp.rem)
    
  }else{
    samp.all <- NA
  }
  
  return(samp.all)
}


# Labels of tips -----------------------------

tips.labels <- function(simpact.trans.net = simpact.trans.net,
                        limitTransmEvents = 7){
  
  seeds.id <- length(simpact.trans.net)
  
  # Add age at sampling
  new.transm.tab <- vector("list", seeds.id)
  
  for(i in 1:seeds.id){
    
    transm.age.i <- as.data.frame(simpact.trans.net[[i]])
    
    age.i <- transm.age.i$SampTime - transm.age.i$TOBRec
    
    transm.age.i <- cbind(transm.age.i, age.i)
    
    new.transm.tab[[i]] <- transm.age.i
    
  }
  
  # ID numbers of Selected networks with at least limitTransmEvents + 1 indiviuals
  
  IDs.transm <- vector()
  
  TransmEventsCountVector <- vector()
  
  for(k in 1:seeds.id){
    trans.net.i.check <- as.data.frame(new.transm.tab[[k]])
    
    if(nrow(trans.net.i.check)>=limitTransmEvents){
      
      TransmEventsCountVector <- c(TransmEventsCountVector, nrow(trans.net.i.check))
      
      IDs.transm <- c(IDs.transm, k)
    }
  }
  
  if(length(IDs.transm)>=1){
    
    ## Binding together all selected transmission transmission networks ##
    
    for (q in 1:length(IDs.transm)){
      
      if(q==1){
        p <- IDs.transm[q]
        trans.sum <- new.transm.tab[[p]]
        rename.id <- paste0(p,".",trans.sum$id,".C")
        trans.sum$id <- rename.id
        trans.sum.rename.id <- trans.sum
      }
      else{
        
        p <- IDs.transm[q]
        
        read.trans.sum <- new.transm.tab[[p]]
        rename.id.read <- paste0(p,".",read.trans.sum$id,".C")
        read.trans.sum$id <- rename.id.read
        trans.sum.rename.id <- rbind(trans.sum.rename.id, read.trans.sum)
      }
      
    }
    
  }
    
    return(trans.sum.rename.id)
    
  }
  