n1 = 60
n2 = 90

# same generation mechanism, except the network is dynamic

  T = 4
  y = array(dim = c(n1 + n2, n1 + n2, T))
  
  latentMean <- matrix(rnorm(6) * 0.5, nrow = 2)
  latentPos1 <- matrix(rnorm(2 * n1) * 0.05 + as.numeric(latentMean), nrow = 2)
  latentPos2 <- matrix(rnorm(2 * n2) * 0.05 + as.numeric(latentMean), nrow = 2)
  
  radius <- matrix(
    c(0.1,0.25,0.25,0.4), nrow = 2
  )
  
  beta <- matrix(
    c(1,0.8,1,1.2), nrow = 2
  )
  for(t in 1:T){
    eta1 <- matrix(0, ncol = n1, nrow = n1)
    for(i in 1:(n1-1)){
      for(j in (i+1):n1){
        HomoDist <- norm(as.matrix(latentPos1[,i] - latentPos1[,j]), 'f')
        eta1[i,j] = beta[1,1] * (radius[1,1] / HomoDist - HomoDist / radius[1,1])
        eta1[j,i] = beta[1,1] * (radius[1,1] / HomoDist - HomoDist / radius[1,1])
      }
    }
    eta1 <- 1/(1 + exp(eta1))
    Y11 <- matrix(runif(n1*n1), nrow = n1)
    Y11 <-  (Y11 > eta1)
    diag(Y11) <- 0
    
    eta2 <- matrix(0, ncol = n2, nrow = n2)
    for(i in 1:(n2-1)){
      for(j in (i+1):n2){
        HomoDist <- norm(as.matrix(latentPos2[,i] - latentPos2[,j]), 'f')
        eta2[i,j] = beta[2,2] * (radius[2,2] / HomoDist - HomoDist / radius[2,2])
        eta2[j,i] = beta[2,2] * (radius[2,2] / HomoDist - HomoDist / radius[2,2])
      }
    }
    eta2 <- 1/(1 + exp(eta2))
    Y22 <- matrix(runif(n2*n2), nrow = n2)
    Y22 <-  (Y22 > eta2)
    diag(Y22) <- 0
    
    
    eta12 <- matrix(0, ncol = n2, nrow = n1)
    eta21 <- matrix(0, ncol = n1, nrow = n2)
    for(i in 1:n1){
      for(j in 1:n2){
        HeteroDist <- norm(as.matrix(latentPos1[,i] - latentPos2[,j]), 'f')
        eta12[i,j] = beta[1,2] * (radius[1,2] / HeteroDist - HeteroDist / radius[1,2])
        eta21[j,i] = beta[2,1] * (radius[2,1] / HeteroDist - HeteroDist / radius[2,1])
      }
    }
    
    eta12 <- 1/(1 + exp(eta12))
    eta21 <- 1/(1 + exp(eta21))
    Y12 <- matrix(runif(n1*n2), nrow = n1)
    Y21 <- matrix(runif(n2*n1), nrow = n2)
    Y12 <-  (Y12 > eta12)
    Y21 <-  (Y21 > eta21)
    
    ## Generated data stroed in Y11, Y12, Y21 in the form of adjacency matricis
    
    y[,,t] <- rbind(cbind(Y11, Y12), cbind(Y21, Y22))
    
    latentPos1 = latentPos1 + matrix(rnorm(2 * n1) * 0.01, nrow = 2)
    latentPos2 = latentPos2 + matrix(rnorm(2 * n2) * 0.01, nrow = 2)
    
  }
  
y1 <- y

# plot true
plot(c(latentPos1[1,],latentPos2[1,]),c(latentPos1[2,],latentPos2[2,]))

setwd("~/Dropbox/Network Data/RPrograms/myModels/DynamicV5")
write.table(y[,,1],"sim_1_t1.csv",quote = F, row.names = F, col.names = F , sep = ',')
write.table(y[,,2],"sim_1_t2.csv",quote = F, row.names = F, col.names = F , sep = ',')
write.table(y[,,3],"sim_1_t3.csv",quote = F, row.names = F, col.names = F , sep = ',')
write.table(y[,,4],"sim_1_t4.csv",quote = F, row.names = F, col.names = F , sep = ',')

node_types = c(rep(0,n1),rep(1,n2))
write.table(node_types,"sim_1_types.csv",quote = F, row.names = F, col.names = F , sep = ',')