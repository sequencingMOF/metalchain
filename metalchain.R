library(dplyr)

rm(list=ls(all=T))


route = 'your path/metalchain/'
knum = 400 # the value of k for k-means algorithm
run_times = 5

#-----------------load data----------------------
setwd(route)
data_original<- read.csv(file="crystal+gaussian noise.csv", header=TRUE, sep="\t")

#-----------------chain searching----------------------
for(run_number in 1:run_times){
 
  data = data_original
  #-----------------create run folder------------------
  setwd(route)
  dir.create(toString(run_number))
  setwd(paste(route,run_number, sep=''))
  
  #-----------------run k means---------------------------------
  cl  = kmeans(data[,1:2], knum, iter.max = 500)
  feature = cbind(data[,-6],cl$cluster) # the 6th column record the k-means cluster a metal belongs to
  data.cluster = feature[order(feature[,3]),] # order from the lowest z to the highest z
  
  #----------------find chain-----------------------------------
  data.cluster = data.cluster[,1:6]
  len = length(data.cluster[,1])
  data.cluster = cbind(data.cluster, matrix(0, nrow=len, ncol=4)) #the 7th to 10th column records
                                                                  #the # of chain in the cluster, 
                                                                  #the # of metal in this chain, 
                                                                  #the number of unenrolled candidates for the chain, 
                                                                  #the number of missing metals for the chain.
  sepa = 0.2253 # the crystallographic vertical distance between two metals in MOF-74
  w = 20 # the anisotropic scalling factor for the penalty of z deviation
  variance.thre = 1 # the threshold of position variance beyond which metals will not be considered for the chain 
  chain.lenthre = 6 # the minimum length of a chain 
  
  for(i in 1:knum){
    chain.ind = which(data.cluster[,6] == i) # get all metals belong to the first k-means cluster
    data.cluster.len = length(chain.ind)
    data.cluster[chain.ind[1],7] = 1 # mark the metal belonging to the first chain of this cluster
    data.cluster[chain.ind[1],8] = 1 # mark the same metal as the first metal of this chain
    zone1.len = 1 # the zone1 length+1; it takes the value of 1 when it is empty.
    zone2.len = 1 # the zone2 length+1; it takes the value of 1 when it is empty.
    variance.lowest = variance.thre
    chain.len = 1 # the current length of the chain
    chain.num = 1 # the # of the chain within the cluster
    candidate.chain = 0 # the number of unenrolled candidate for a chain
    k=1 # k is the index of the last enrolled point of the chain
    j=2 # j is the index of the candidate under evaluation
    xc = data.cluster[chain.ind[k],1] #initialize the center of the chain
    yc = data.cluster[chain.ind[k],2] #initialize the center of the chain
    zc = data.cluster[chain.ind[k],3] #initialize the center of the chain
    miss = 0 # means the number of missing atoms in the chain
    while(j < (data.cluster.len + 1)){
      if((data.cluster[chain.ind[j],3]-zc)<1.5*sepa){    #judge if the metal under evaluation is in zone #1  
        j=j+1 #prepare to evaluate the next metal
        temp = (data.cluster[chain.ind[j-1],1]-xc)^2 + (data.cluster[chain.ind[j-1],2]-yc)^2 + w*(data.cluster[chain.ind[j-1],3]-zc-sepa)^2 # the position variance
        if(temp < variance.lowest){ #judge if the metal under evaluation is by far the most fit candidate
          variance.lowest = temp 
          zone.fittest = j-1
          zone1.len = zone1.len + 1
        }
      }
      else{ #the metal is not in zone #1
        if(zone1.len == 1){ # the metal is not in zone #1 and zone #1 is empty
          if(abs(data.cluster[chain.ind[j],3]-zc-2*sepa)<0.5*sepa){ #judge if the metal under evaluation is in zone #2  
            j=j+1 #prepare to evaluate the next metal
            temp = (data.cluster[chain.ind[j-1],1]-xc)^2 + (data.cluster[chain.ind[j-1],2]-yc)^2 + w*(data.cluster[chain.ind[j-1],3]-zc-2*sepa)^2
            if(temp < variance.lowest){ #judge if the metal under evaluation is by far the most fit candidate
              variance.lowest = temp
              zone.fittest = j-1
              zone2.len = zone2.len + 1
            }
          }
          else{ # the metal is neither in zone #1 or zone #2, 
            if(zone2.len == 1){ # both zone #1 and zone #2 are empty, terminate this chain
              if(data.cluster[chain.ind[k],8] < chain.lenthre){ # the chain is too short
                for (t in (k-candidate.chain-chain.len+1):(k)) { # the premature chain is removed
                  data.cluster[chain.ind[t],7] = 0
                  data.cluster[chain.ind[t],8] = 0
                }
              }
              else{ # the chain is long enough and can be harvested
                chain.num = chain.num+1 
                for (t in (k-candidate.chain-chain.len+1):(k)) {
                  data.cluster[chain.ind[t],9] = candidate.chain
                  data.cluster[chain.ind[t],10] = miss
                }
              }
              k=k+1 #initialize the search for the next chain
              j=k+1 
              data.cluster[chain.ind[k],7] = chain.num
              data.cluster[chain.ind[k],8] = 1
              xc = data.cluster[chain.ind[k],1]
              yc = data.cluster[chain.ind[k],2]
              zc = data.cluster[chain.ind[k],3]
              chain.len = 1
              variance.lowest = variance.thre
              candidate.chain = 0
              zone1.len = 1
              zone2.len = 1
              miss = 0
            }
            else{ # zone #1 is empty but zone #2 was filled
              data.cluster[chain.ind[zone.fittest],7] = data.cluster[chain.ind[k],7] # enroll the best candidate in zone #2 
              data.cluster[chain.ind[zone.fittest],8] = data.cluster[chain.ind[k],8] + 2
              xc = (xc*chain.len+data.cluster[chain.ind[zone.fittest],1])/(chain.len+1) # update the center of the chain
              yc = (yc*chain.len+data.cluster[chain.ind[zone.fittest],2])/(chain.len+1)
              zc = (data.cluster[chain.ind[zone.fittest],3]-zc-2*sepa)/(chain.len+1) + zc + 2*sepa
              chain.len = chain.len + 1
              candidate.chain = candidate.chain + zone.fittest - k - 1 #update the number of candidate in the chain
              k = zone.fittest #continue evaluation from the next metal after the enrolled 
              zone1.len = 1 #initialize the search for the next metal
              zone2.len = 1
              variance.lowest = variance.thre
              j=k+1
              miss = miss + 1
            }
          }
        }
        else{ # the metal is not in zone #1 and zone #1 is filled
          data.cluster[chain.ind[zone.fittest],7] = data.cluster[chain.ind[k],7] # enroll the best candidate in zone #1 
          data.cluster[chain.ind[zone.fittest],8] = data.cluster[chain.ind[k],8] + 1
          xc = (xc*chain.len+data.cluster[chain.ind[zone.fittest],1])/(chain.len+1)
          yc = (yc*chain.len+data.cluster[chain.ind[zone.fittest],2])/(chain.len+1)
          zc = (data.cluster[chain.ind[zone.fittest],3]-zc-sepa)/(chain.len+1) + zc + sepa
          chain.len = chain.len + 1
          candidate.chain = candidate.chain + zone.fittest - k - 1
          k = zone.fittest
          zone1.len = 1
          zone2.len = 1
          variance.lowest = variance.thre
          j=k+1
        }
      }
    }
    if(data.cluster[chain.ind[k],8] < chain.lenthre){ #finish the harvest of the last chain
      for (t in (k-candidate.chain-chain.len+1):(k)) {
        data.cluster[chain.ind[t],7] = 0
        data.cluster[chain.ind[t],8] = 0
      }
    }
    else{
      for (t in (k-candidate.chain-chain.len+1):(k)) {
        data.cluster[chain.ind[t],9] = candidate.chain
        data.cluster[chain.ind[t],10] = miss
      }
    }
  }
  
  raw.chain = data.cluster[order(data.cluster[,6]),] # order by  cluster
  raw.chain = raw.chain[which(raw.chain[,7]!=0),] #remove unassigned metals
  
   #--------------calculate within-chain position variance-----------------------------
  
  raw.chain = raw.chain[,1:10]
  raw.chain.len = length(raw.chain[,1])
  raw.chain = cbind(raw.chain, matrix(0, nrow=raw.chain.len, ncol=1)) #add a column to record chain variance
  
  chain.len = 1
  xc = raw.chain[1,1]
  yc = raw.chain[1,2]
  zc = raw.chain[1,3]
  chain.value = 1 
  cluster.value = raw.chain[1,6]
  
  for(i in 1:(raw.chain.len-1)){
    if(raw.chain[i+1,7] == chain.value && raw.chain[i+1,6] == cluster.value){
      chain.len = chain.len + 1
      xc = xc + raw.chain[i+1,1]
      yc = yc + raw.chain[i+1,2]
      zc = zc + raw.chain[i+1,3] - sepa * (raw.chain[i+1,8] - 1) 
    }
    else{
      xc = xc/chain.len
      yc = yc/chain.len
      zc = zc/chain.len
      chain.vari = 0
      for(j in (i-chain.len+1):i){
        chain.vari = chain.vari + (raw.chain[j,1]-xc) **2 + (raw.chain[j,2]-yc) **2 + (raw.chain[j,3]-zc-(raw.chain[j,8]-1) * sepa) **2 *w
      }
      chain.vari = chain.vari/chain.len
      for(j in (i-chain.len+1):i){
        raw.chain[j,11] = chain.vari # record the within-chain position variance 
      }
      xc = raw.chain[i+1,1] #initialize the measuring the next chain
      yc = raw.chain[i+1,2]
      zc = raw.chain[i+1,3]
      chain.len = 1
      chain.value = raw.chain[i+1,7]
      cluster.value = raw.chain[i+1,6]
    }
  }
  
  xc = xc/chain.len
  yc = yc/chain.len
  zc = zc/chain.len
  chain.vari = 0
  for(j in (i-chain.len+2):(i+1)){
    chain.vari = chain.vari + (raw.chain[j,1]-xc)^2 + (raw.chain[j,2]-yc)^2 + (raw.chain[j,3] - zc - (raw.chain[j,8]-1) * sepa)^2 * w
  }
  chain.vari = chain.vari/chain.len
  
  for(j in (i-chain.len+2):(i+1)){
    raw.chain[j,11] = chain.vari
  }
  
  #-------------------------Selection of chains------------------------------------------------------
  dat = data.frame(raw.chain)
  colnames(dat) = c("x","y","z","chain_crys","atom_crys","cluster","chain.num","atom.num","candidate","miss","variance")
  chain.summary = summarise(group_by(dat, cluster, chain.num, candidate, miss, variance), np= n())
  
  chain.summary[is.na(chain.summary)] = 0
  
  var.thre = 0.15 # the threshold of witin-chain variance
  atom.thre = 0.5 #the threshold of the confidence index
  chain.summary$confidence <- (1-(chain.summary$candidate + chain.summary$miss) / chain.summary$np) # the confidence index considers 
                                                                                                    # the number of unenrolled candidates
                                                                                                    # and missing atoms
  
  ind = which(chain.summary$variance < var.thre & (chain.summary$confidence > atom.thre) )
  good.chain = chain.summary[ind,] #chains that passed the selection are recorded in good.chain
  
  #----------------------get all the metals belonging to the selected chains---------------------------------------
  chain.num = 0
  s = 1
  k = 1
  while(k < length(raw.chain[,1])+1 && s < length(good.chain$cluster)+1){
    if(raw.chain[k,6]==good.chain[s,1]&&raw.chain[k,7]==good.chain[s,2]){
      chain.num = rbind(chain.num, k)
      k = k+1
    }
    else{
      if(s == length(good.chain$cluster)){break}
      if(raw.chain[k,6]==good.chain[s+1,1]&&raw.chain[k,7]==good.chain[s+1,2]){s=s+1}
      else{k=k+1}
    }
  }
  chain.num = chain.num[-1]
  Final.Chain = raw.chain[chain.num,1:8]
  Final.Chain = as.data.frame(Final.Chain)
  raw.Final.Chain = raw.chain[,1:8]
  raw.Final.Chain = as.data.frame(raw.Final.Chain)
  
  colnames(Final.Chain) = c("x","y","z","rod#","position on rod","cluster#","chain#","position on chain")
  colnames(raw.Final.Chain) = c("x","y","z","rod#","position on rod","cluster#","chain#","position on chain")

  write.table(Final.Chain, file="atom info.csv",row.names = FALSE, sep = "\t") # the metals belonging to the selected chains
  write.table(raw.Final.Chain, file="atom info_all.csv",row.names = FALSE,sep = "\t") # all the metals belonging to all chains
  
  #----------------------------------------------------------------------------------------------------
  #-------------------analyze error of selected chains-------------------------
  cluster.value = Final.Chain[1,6]
  chain.value = Final.Chain[1,7]
  original.chain = matrix(0, ncol=2, nrow=1)
  original.chain[1,1] = Final.Chain[1,4] #the first rod in the crystal
  original.chain[1,2] = 1 #the second column records the number of metal has been found for this rod
  error.stat = matrix(0, ncol=4,nrow=0) #record error
  
  for(i in 1:(length(Final.Chain[,1])-1)){
    if(Final.Chain[i+1,7] == chain.value && Final.Chain[i+1,6] == cluster.value){ #judge if the metal under evaluation is assigned to the chain
      found = 0
      for (j in 1:length(original.chain[,1])){
        if(Final.Chain[i+1,4] == original.chain[j,1]){ # judge if the metal under evaluation is correctly assigned 
          original.chain[j,2] = original.chain[j,2] + 1
          found = 1
        }
      }    
      if (found == 0){ # the metal under evaluation belongs to a rod that has not been evalauted yet
        new.original = matrix(0, ncol = 2, nrow=1)
        new.original[1,1] = Final.Chain[i+1,4] #evaluat this rod
        new.original[1,2] = 1
        original.chain = rbind(original.chain, new.original)
      }
    }
    else{ # the metal was not assigned to the current chain and was assigned as the first metal of the next chain
      new.error = matrix(0, ncol=4,nrow=1)
      new.error[1,1] = cluster.value
      new.error[1,2] = chain.value
      new.error[1,3] = apply(original.chain,2,max)[2]
      new.error[1,4] = colSums(original.chain)[2]
      error.stat = rbind(error.stat, new.error)
      chain.value = Final.Chain[i+1,7]
      cluster.value = Final.Chain[i+1,6]
      original.chain = matrix(0, ncol=2, nrow=1)
      original.chain[1,1] = Final.Chain[i+1,4]
      original.chain[1,2] = 1
    }
  }
  
  
  new.error = matrix(0, ncol=4,nrow=1)
  new.error[1,1] = cluster.value
  new.error[1,2] = chain.value
  new.error[1,3] = apply(original.chain,2,max)[2]
  new.error[1,4] = colSums(original.chain)[2]
  error.stat = rbind(error.stat, new.error)
  
  error = data.frame(error.stat)
  colnames(error) = c("cluster#","chain#","correct","total")
  write.table(error, file="error.csv",row.names = FALSE, sep = "\t")
  #----------------------------------------------------------------------
  
  
  #-------------------analyze error of all chains-------------------------
  cluster.value = raw.Final.Chain[1,6]
  chain.value = raw.Final.Chain[1,7]
  original.chain = matrix(0, ncol=2, nrow=1)
  original.chain[1,1] = raw.Final.Chain[1,4]
  original.chain[1,2] = 1
  error.stat = matrix(0, ncol=4,nrow=0)
  
  for(i in 1:(length(raw.Final.Chain[,1])-1)){
    if(raw.Final.Chain[i+1,7] == chain.value && raw.Final.Chain[i+1,6] == cluster.value){
      found = 0
      for (j in 1:length(original.chain[,1])){
        if(raw.Final.Chain[i+1,4] == original.chain[j,1]){
          original.chain[j,2] = original.chain[j,2] + 1
          found = 1
        }
      }    
      if (found == 0){
        new.original = matrix(0, ncol = 2, nrow=1)
        new.original[1,1] = raw.Final.Chain[i+1,4]
        new.original[1,2] = 1
        original.chain = rbind(original.chain, new.original)
      }
    }
    else{
      new.error = matrix(0, ncol=4,nrow=1)
      new.error[1,1] = cluster.value
      new.error[1,2] = chain.value
      new.error[1,3] = apply(original.chain,2,max)[2]
      new.error[1,4] = colSums(original.chain)[2]
      error.stat = rbind(error.stat, new.error)
      chain.value = raw.Final.Chain[i+1,7]
      cluster.value = raw.Final.Chain[i+1,6]
      original.chain = matrix(0, ncol=2, nrow=1)
      original.chain[1,1] = raw.Final.Chain[i+1,4]
      original.chain[1,2] = 1
    }
  }
  
  
  new.error = matrix(0, ncol=4,nrow=1)
  new.error[1,1] = cluster.value
  new.error[1,2] = chain.value
  new.error[1,3] = apply(original.chain,2,max)[2]
  new.error[1,4] = colSums(original.chain)[2]
  error.stat = rbind(error.stat, new.error)
  
  error = data.frame(error.stat)
  colnames(error) = c("cluster#","chain#","correct","total")
  write.table(error, file="error_all.csv",row.names = FALSE, sep = "\t")
  #----------------------------------------------------------------------
  
}

setwd(route)


#--------------------------------------------------------summary error for all the runs ---------------------------
error_total = matrix(0, nrow=run_times, ncol = 3)

for(i in 1:run_times){
  setwd(paste(route,i, sep=''))
  error <- read.csv(file="error.csv", header=TRUE, sep="\t")
  error_total[i,1] = colSums(error)[3]
  error_total[i,2] = colSums(error)[4]
} 
error_total[,3] = error_total[,1]/error_total[,2]
mean = colMeans(error_total)
error_total = rbind(error_total, mean)
setwd(route)
colnames(error_total) = c("correct","total","accuracy")
write.table(error_total, file="error_stat.csv",row.names = FALSE,sep = "\t")
#-----------------------------------------------------------------------------------------------------

#--------------------------------------------------------summary of error_all for all the runs ---------------------------
error_total = matrix(0, nrow=run_times, ncol = 3)

for(i in 1:run_times){
  setwd(paste(route,i, sep=''))
  error <- read.csv(file="error_all.csv", header=TRUE, sep="\t")
  error_total[i,1] = colSums(error)[3]
  error_total[i,2] = colSums(error)[4]
} 
error_total[,3] = error_total[,1]/error_total[,2]
mean = colMeans(error_total)
error_total = rbind(error_total, mean)
setwd(route)
colnames(error_total) = c("correct","total","accuracy")
write.table(error_total, file="error_all_stat.csv",row.names = FALSE,sep = "\t")
#---------------------------------------------------------------------------------------------------------  