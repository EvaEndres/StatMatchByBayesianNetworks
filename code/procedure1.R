# fit Bayesian networks

# read data
source("read_data.R")

# load example data
load("../data/allbus_dat.RData")

# load packages
library(bnlearn)
library(colorspace)
library(gtools)
library(tikzDevice)

# load code
source("functions.R")

# names of variables
common <- names(datA)[names(datA) %in% names(datB)]
specificA <- names(datA)[!(names(datA) %in% names(datB))]
specificB <- names(datB)[!(names(datB) %in% names(datA))]


######################################
# step 1, procedure 1: fit DAG for X #
######################################

# summarize observations of the common variables in X
X <- rbind(datA[,common], datB[,common])
dim(X)

# learn the DAG
X.bn.analysis <- bn.analysis(data=X, black=NULL, white=NULL, modelstring=NULL, col.code=rep(rainbow_hcl(3)[1], times=ncol(X)), h.nodes=sort(common))

# save blacklist for X (edges that must not appear in the final DAG)
blackX <- black(bn.an=X.bn.analysis, data=X)

# save structure of common variables in whitelist (egdes that must appear in the final DAG)
whiteX <- X.bn.analysis$avg.boot.hc$arcs
whiteX


##############################################
# step 1, procedure 1: fit DAG on data set A #
##############################################

# ensure fork connection and pass X-structure as prior knowledge
blackA <- matrix(c(rep(specificA, each=length(common)), rep(common, length(specificA))), ncol=2)
blackA <- rbind(blackA, blackX)

# learn the DAG
A.bn.analysis <- bn.analysis(data=datA, black=blackA, white=whiteX, col.code=c(rep(rainbow_hcl(3)[1], times=length(common)), rep(rainbow_hcl(3)[2], times=length(specificA))), modelstring=NULL, h.nodes=c(common, specificA))
warnings() # cylcles occured during model averaging!

# fit DAG again without cycle-causing arc
emptyA <- empty.graph(c(common, specificA))
arcs(emptyA) <- arcs(A.bn.analysis$avg.boot.hc)
msA <- modelstring(emptyA)
A.bn.analysis <- bn.analysis(data=datA, black=NULL, white=NULL, col.code=c(rep(rainbow_hcl(3)[1], times=length(common)), rep(rainbow_hcl(3)[2], times=length(specificA))), modelstring=msA, h.nodes=c(common, specificA))

# save structure of A in a whitelist
whiteA <- A.bn.analysis$avg.boot.hc$arcs
whiteA

##############################################
# step 1, procedure 1: fit DAG on data set B #
##############################################

# ensure fork connection and pass X-structure as prior knowledge
blackB <- matrix(c(rep(specificB, each=length(common)), rep(common, length(specificB))), ncol=2)
blackB <- rbind(blackB, blackX)

# learn the DAG
B.bn.analysis <- bn.analysis(data=datB, black=blackB, white=whiteX, col.code=c(rep(rainbow_hcl(3)[1], times=length(common)), rep(rainbow_hcl(3)[3], times=length(specificB))), modelstring=NULL, h.nodes=c(common, specificB))
warnings()

# cylcles occured during model averaging!
# fit DAG again without cycle-causing arc
# emptyB <- empty.graph(c(common, specificB))
# arcs(emptyB) <- arcs(B.bn.analysis$avg.boot.hc)
# msB <- modelstring(emptyB)
# B.bn.analysis <- bn.analysis(data=datB, black=NULL, white=NULL, col.code=c(rep(rainbow_hcl(3)[1], times=length(common)), rep(rainbow_hcl(3)[3], times=length(specificB))), modelstring=msB, h.nodes=c(common, specificB))


# save structure of specificB variables in whitelist
whiteB <- B.bn.analysis$avg.boot.hc$arcs
whiteB

#############################################
# step 1, procedure 1: merge DAG structures #
#############################################

bnXYZ <- empty.graph(c(common, specificA, specificB))
arcs(bnXYZ) <- rbind(whiteX, whiteA, whiteB)[which(!duplicated(rbind(whiteX, whiteA, whiteB))),]

# plot DAG with the Graphviz package
h.nodes = c(common, specificA,  specificB)
col.code <- c(rep(rainbow_hcl(3)[1], times=length(common)), rep(rainbow_hcl(3)[2], times=length(specificA)), rep(rainbow_hcl(3)[3], times=length(specificB)))

high = list(nodes = h.nodes,
            col = col.code, 
            fill = col.code, 
            textCol = c(rep("black", times=length(common)), rep("black", times=(length(specificA)+length(specificB)))))

graphviz.plot(bnXYZ, shape = "ellipse", highlight = high)

###############################
# step 2: estimate parameters #
###############################

dist <- list()
length(dist) <- length(c(common, specificA, specificB))
names(dist) <- c(common, specificA, specificB)

for(j in 1:length(dist)){
  
  if(names(dist)[j] %in% common){bn.fit <- X.bn.analysis$fit}
  if(names(dist)[j] %in% specificA){bn.fit <- A.bn.analysis$fit}
  if(names(dist)[j] %in% specificB){bn.fit <- B.bn.analysis$fit}
  dist[[j]] <- bn.fit[[which(names(bn.fit)==names(dist)[j])]]$prob
  if(any(is.na(dist[[j]]))){
    dist[[j]] <- NA.unif(dist[[j]])
  }
  
}
match.fit <- custom.fit(bnXYZ, dist = dist)

save(match.fit, file="../data/match.fit.rda")

#######################################
# step 3: generate synthetic data set #
#######################################

# impute Z in A
matchA <- data.frame(datA, doctor=NA, hospital=NA, smoke=NA, alcohol=NA, health.state=NA)

# topological ordering: doctor, hospital, alcohol, smoke, health.state
for(i in 1:nrow(matchA)){
  set.seed(10+(100*i))
  matchA$doctor[i] <- sample(levels(datB$doctor),1,prob=match.fit$doctor$prob[,which(levels(datB$sex)==matchA$sex[i])])
  matchA$hospital[i] <- sample(levels(datB$hospital),1,prob=match.fit$hospital$prob[,which(levels(datB$sex)==matchA$sex[i])
                                                                                    ,which(levels(datB$employed)==matchA$employed[i])
                                                                                    ,which(levels(datB$marital.status)==matchA$marital.status[i])
                                                                                    ,which(levels(datB$income)==matchA$income[i])
                                                                                    ,which(levels(datB$doctor)==matchA$doctor[i])
                                                                                    ])
  matchA$alcohol[i] <- sample(levels(datB$alcohol),1,prob=match.fit$alcohol$prob[,which(levels(datB$sex)==matchA$sex[i])
                                                                                 ,which(levels(datB$marital.status)==matchA$marital.status[i])
                                                                                 ,which(levels(datB$hospital)==matchA$hospital[i])
                                                                                 ])
  matchA$smoke[i] <- sample(levels(datB$smoke),1,prob=match.fit$smoke$prob[,which(levels(datB$sex)==matchA$sex[i])
                                                                           ,which(levels(datB$hospital)==matchA$hospital[i])
                                                                           ])
  matchA$health.state[i] <- sample(levels(datB$health.state),1,prob=match.fit$health.state$prob[,which(levels(datB$sex)==matchA$sex[i])
                                                                                                ,which(levels(datB$doctor)==matchA$doctor[i])
                                                                                                ,which(levels(datB$hospital)==matchA$hospital[i])
                                                                                                ])
}

# impute Y in B
matchB <- data.frame(datB, confession=NA, church=NA, exp.god=NA, exp.supernat=NA, faith.healer=NA, pray=NA)

# topological ordering: pray, confession, church, exp.god, exp.supernat, faith.healer
for(i in 1:nrow(matchB)){
  set.seed(10+(100*i))
  matchB$pray[i] <- sample(levels(datA$pray),1,prob=match.fit$pray$prob)
  matchB$confession[i] <- sample(levels(datA$confession),1,prob=match.fit$confession$prob[,which(levels(datA$sex)==matchB$sex[i])
                                                                                          ,which(levels(datA$graduation)==matchB$graduation[i])
                                                                                          ,which(levels(datA$pray)==matchB$pray[i])
                                                                                          ])
  matchB$church[i] <- sample(levels(datA$church),1,prob=match.fit$church$prob[,which(levels(datA$confession)==matchB$confession[i])])
  matchB$exp.god[i] <- sample(levels(datA$exp.god),1,prob=match.fit$exp.god$prob[,which(levels(datA$sex)==matchB$sex[i])
                                                                                 ,which(levels(datA$confession)==matchB$confession[i])
                                                                                 ,which(levels(datA$church)==matchB$church[i])
                                                                                 ,which(levels(datA$pray)==matchB$pray[i])
                                                                                 ])
  matchB$exp.supernat[i] <- sample(levels(datA$exp.supernat),1,prob=match.fit$exp.supernat$prob[,which(levels(datA$church)==matchB$church[i])
                                                                                                ,which(levels(datA$exp.god)==matchB$exp.god[i])
                                                                                                ])
  matchB$faith.healer[i] <- sample(levels(datA$faith.healer),1,prob=match.fit$faith.healer$prob[,which(levels(datA$sex)==matchB$sex[i])
                                                                                                ,which(levels(datA$graduation)==matchB$graduation[i])
                                                                                                ,which(levels(datA$marital.status)==matchB$marital.status[i])
                                                                                                ,which(levels(datA$church)==matchB$church[i])
                                                                                                ,which(levels(datA$exp.god)==matchB$exp.god[i])
                                                                                                ,which(levels(datA$exp.supernat)==matchB$exp.supernat[i])
                                                                                                ])
}

# combine imputed A and B
match.dat_post <- rbind(matchA[,c(common, specificA, specificB)], matchB[,c(common, specificA, specificB)])
match.dat_post$doctor <- factor(match.dat_post$doctor, levels=levels(datB$doctor))
match.dat_post$hospital <- factor(match.dat_post$hospital, levels=levels(datB$hospital))
match.dat_post$smoke <- factor(match.dat_post$smoke, levels=levels(datB$smoke))
match.dat_post$alcohol <- factor(match.dat_post$alcohol, levels=levels(datB$alcohol))
match.dat_post$health.state <- factor(match.dat_post$health.state, levels=levels(datB$health.state))


save(match.dat_post, file="../data/match_dat_post.rda")

