######################################
# comparison with original GGSS data #
######################################

rm(list=ls())

# load synthetic data and original data
load("../data/match_dat_post.rda")
load("../data/match.dat.union_post.rda")
load("../data/match.dat.intersection_post.rda")
load("../data/allbus.compl.rda")
allbus <- allbus.compl
match <- match.dat_post

# synthetic data sets are equal for all introduced matching procedures
all(match.dat_post==match.dat.union_post)
all(match.dat_post==match.dat.intersection_post)

# load source code
source("functions.R")

# names of variables
load("../data/allbus_dat.RData")
common <- names(datA)[names(datA) %in% names(datB)]
specificA <- names(datA)[!(names(datA) %in% names(datB))]
specificB <- names(datB)[!(names(datB) %in% names(datA))]

###############################
# comparison of the marginals #
###############################

spec <- list()
length(spec) <- length(c(specificA, specificB))
names(spec) <- c(specificA, specificB)
res.jsd <- res.p <- spec

for(j in 1:length(spec)){
  
  res.jsd[[j]] <- jsd(table(allbus[which(names(allbus)==names(res.jsd)[j])])/nrow(allbus), table(match[which(names(match)==names(res.jsd)[j])])/nrow(match))
  
  test.data <- data.frame(
    rbind(match[which(names(match)==names(res.p)[j])], allbus[which(names(allbus)==names(res.p)[j])]), 
    id=c(rep("match", nrow(match)), rep("allbus", nrow(allbus)))
  )
  test.data[,1] <- as.factor(test.data[,1])
  res.p[[j]] <- chisq.test(test.data[,1], test.data[,2])$p.value
}

# results of the Jensen-Shannon divergence
save(res.jsd, file="../data/res.jsd.rda")
round(unlist(res.jsd),3)

# p values of the univariate chi^2 test
save(res.p, file="../data/res.p.rda")
round(unlist(res.p),3)

##############################################
# comparison of the contingency coefficients #
##############################################

library(DescTools)

#matrix of bivariate contingency coefficients of all specific variables in the original GGSS data
cont.matrix.allbus <- matrix(0, ncol=(length(specificB)+length(specificA)), nrow=(length(specificB)+length(specificA)))
colnames(cont.matrix.allbus) <- rownames(cont.matrix.allbus) <- c(specificA, specificB)

for(i in 1:nrow(cont.matrix.allbus)){
  for(j in 1:ncol(cont.matrix.allbus)){
    cont.matrix.allbus[i,j] <- ContCoef(table(allbus[,which(names(allbus)==colnames(cont.matrix.allbus)[j])], 
                                              allbus[,which(names(allbus)==rownames(cont.matrix.allbus)[i])]), correct = TRUE)
  }
}


#matrix of bivariate contingency coefficients of all specific variables in the original GGSS data
cont.matrix.match <- matrix(0, ncol=(length(specificB)+length(specificA)), nrow=(length(specificB)+length(specificA)))
colnames(cont.matrix.match) <- rownames(cont.matrix.match) <- c(specificA, specificB)

for(i in 1:nrow(cont.matrix.match)){
  for(j in 1:ncol(cont.matrix.match)){
    cont.matrix.match[i,j] <- ContCoef(table(match[,which(names(match)==colnames(cont.matrix.match)[j])], 
                                             match[,which(names(match)==rownames(cont.matrix.match)[i])]), correct = TRUE)
  }
}

#compute absolute deviations
dev.match.post <- abs(cont.matrix.allbus-cont.matrix.match)
dev.match.post[lower.tri(dev.match.post)] <- NA
diag(dev.match.post) <- NA
dev.match.post
mean(dev.match.post, na.rm=TRUE)
min(dev.match.post, na.rm=TRUE)
max(dev.match.post, na.rm=TRUE)

mean(dev.match.post, na.rm=TRUE)
sd(dev.match.post, na.rm=TRUE)

