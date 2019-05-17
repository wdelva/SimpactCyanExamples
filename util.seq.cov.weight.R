# Function to tranform dates in named vector to be handled by treedater
# calendar.dates File containing named sampling calendar times for each sequence
# count.start Calendar year when the simulation started
# endsim Number of years when the simulation was done
dates.Transform.NamedVector  <- function(dates=dates,
                                         count.start = 1977,
                                         endsim = 40){
  
  dates.val <- endsim - dates$V2 + count.start # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  
  return(dates.val)
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