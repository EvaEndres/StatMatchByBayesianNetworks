# perform structure learning and plot the resulting structure

bn.analysis <- function(data, black=NULL, white=NULL, col.code=NULL, h.nodes=NULL, modelstring=NULL){
  
  if(is.null(modelstring)){
    #learn 500 network structures on the common variables from bootstrap samples 
    set.seed(12345)
    boot.hc <- boot.strength(data=data, m=round(2*nrow(data)/3), R=500, algorithm="hc", algorithm.args=list(blacklist=black, whitelist=white), cpdag = FALSE)
    
    #average the bootstrap-networks and set threshold=0
    avg.boot.hc <- averaged.network(boot.hc, threshold=0)
  }
  else{
    avg.boot.hc <- model2network(modelstring)
  }
  
  #plot with the Graphviz package
  plot.net(avg.boot.hc=avg.boot.hc, col.code=col.code, h.nodes=h.nodes)
  
  #estimate parameters
  fit <- bn.fit(avg.boot.hc, data)
  
  #save structure of common variables in whitelist
  whiteX <- avg.boot.hc$arcs
  
  #return
  if(is.null(modelstring)){
    return(list(boot.hc=boot.hc, avg.boot.hc=avg.boot.hc, fit=fit))
  }
  else{
    return(list(avg.boot.hc=avg.boot.hc, fit=fit))
  }
}


# generate new blacklist for data

black <- function(bn.an, data){
  whiteX <- bn.an$avg.boot.hc$arcs
  bmat <- combinations(n = ncol(data), r = 2, v = common, repeats.allowed = FALSE)
  bmat <- rbind(bmat, bmat[,c(2,1)])
  bmat <- bmat[order(bmat[,1]),]
  del <- vector()
  for(i in 1:nrow(whiteX)){
    del <- c(del,which((bmat[,1]==whiteX[i,1])&((bmat[,2]==whiteX[i,2]))))
  }
  blackX <- bmat[-del,]
  return(blackX)
}

plot.net <- function(avg.boot.hc, col.code=NULL, h.nodes=NULL){
  
  #plot with the Graphviz package
  high = list(nodes = h.nodes, col = col.code,  fill = col.code, textCol = "black") 
  graphviz.plot(avg.boot.hc, shape = "ellipse", highlight = high)
  
}

# assume uniform distribution for child-nodes given rare parent instantiations

NA.unif <- function(cpt){
  NA.index <- which(is.na(cpt), arr.ind=TRUE)
  n.levels <- vector(length=length(dimnames(cpt)))
  
  for(j in 1:length(n.levels)){
    n.levels[j] <- length(dimnames(cpt)[[j]])
  }
  
  n.run <- nrow(NA.index)/n.levels[[1]]
  
  NA.run <- which(is.na(cpt), arr.ind=FALSE)
  cpt[NA.run] <- 1/n.levels[1]
  
  return(cpt)
  
}

# function to compute the Jensen-Shannon divergence

jsd <- function(p,q){
  
  require(entropy)
  
  m <- 0.5*(p+q)
  pm <- KL.empirical(p, m, unit="log2")
  qm <- KL.empirical(p, m, unit="log2")
  jsd <- 0.5*pm+0.5*qm
  return(jsd)
  
  
}

# edge combination

combi <- function(whiteA, whiteB, method=c("union","intersection")){
  
  idA <- which((whiteA[,1]%in% common)&(whiteA[,2]%in% common))
  idB <- which((whiteB[,1]%in% common)&(whiteB[,2]%in% common))
  
  A <- whiteA[idA,]
  B <- whiteB[idB,]
  
  specA <- whiteA[-idA,]
  specB <- whiteB[-idB,]
  
  mat <- rbind(A, B)
  dup <- mat[duplicated(mat),]
  arcs <- mat[!duplicated(mat),]
  
  if(method=="intersection"){
    return(list(all.arcs=rbind(dup, specA, specB), arcs.common=dup))
  }
  if(method=="union"){
    return(list(all.arcs=rbind(arcs, specA, specB)))
  }
  
}