# The GGSS data can be downloaded from https://dbk.gesis.org/DBKSearch/SDESC2.asp?no=4614&tab=3&db=D&dab=0

# read data
library(foreign)
options(warn=-1) # http://stackoverflow.com/questions/21228518/what-is-attr-value-labels-when-reading-spss-into-r
allbus <-  read.spss("../data/ZA4614_v1-1-1.sav",  to.data.frame=TRUE, use.value.labels=TRUE)
options(warn=0)
 
# extract demographic variables
demo <- allbus[,c("V217","V220","V230","V622","V274","V344")]
names(demo) <- c("sex", "age", "graduation", "employed", "marital.status", "income")

# extract religious variables
rel <- allbus[,c("V188","V269","V143","V144","V163","V271")]
names(rel) <- c("confession","church","exp.god","exp.supernat","faith.healer","pray")

# extract health-related variables
health <- allbus[,c("V574","V576","V587","V588","V591")]
names(health) <- c("doctor","hospital","smoke","alcohol","health.state")

# complete data
allbus.compl <- na.omit(cbind(demo, rel, health))

# discretize continuous variables
library(bnlearn)
?discretize
allbus.discrete <- discretize(allbus.compl, method="interval", breaks=20)

# exclude categories with less than 20 observations
allbus.discrete$age[allbus.discrete$age %in% levels(allbus.discrete$age)[summary(allbus.discrete$age)<20]] <- NA
allbus.discrete$graduation[allbus.discrete$graduation%in% levels(allbus.discrete$graduation)[summary(allbus.discrete$graduation)<20]] <- NA
allbus.discrete$employed[allbus.discrete$employed%in% levels(allbus.discrete$employed)[summary(allbus.discrete$employed)<20]] <- NA
allbus.discrete$income[allbus.discrete$income%in% levels(allbus.discrete$income)[summary(allbus.discrete$income)<20]] <- NA
allbus.discrete$marital.status[allbus.discrete$marital.status%in% levels(allbus.discrete$marital.status)[summary(allbus.discrete$marital.status)<20]] <- NA

allbus.discrete$church[allbus.discrete$church%in% levels(allbus.discrete$church)[summary(allbus.discrete$church)<20]] <- NA
allbus.discrete$confession[allbus.discrete$confession%in% levels(allbus.discrete$confession)[summary(allbus.discrete$confession)<20]] <- NA

allbus.discrete$smoke[allbus.discrete$smoke%in% levels(allbus.discrete$smoke)[summary(allbus.discrete$smoke)<20]] <- NA
allbus.discrete$alcohol[allbus.discrete$alcohol%in% levels(allbus.discrete$alcohol)[summary(allbus.discrete$alcohol)<20]] <- NA


# complete data
allbus.compl <- na.omit(allbus.discrete)

# relevel some variables with zero observations in some categories
allbus.compl$graduation <- factor(allbus.compl$graduation)
allbus.compl$employed <- factor(allbus.compl$employed)
allbus.compl$marital.status <- factor(allbus.compl$marital.status)
allbus.compl$age <- factor(allbus.compl$age)
allbus.compl$income <- factor(allbus.compl$income)

allbus.compl$church <- factor(allbus.compl$church)
allbus.compl$confession <- factor(allbus.compl$confession)

allbus.compl$smoke <- factor(allbus.compl$smoke)
allbus.compl$alcohol <- factor(allbus.compl$alcohol)


# balanced allocation in two data set
n.A <- round(nrow(allbus.compl)/2)
set.seed(12345)
index.A <- sample(1:nrow(allbus.compl),n.A, replace=FALSE)
datA <- allbus.compl[index.A,c(names(demo), names(rel))]
datB <- allbus.compl[-index.A,c(names(demo), names(health))]

#save data
save(datA, datB, file="../data/allbus_dat.RData")
save(allbus.compl, file="../data/allbus.compl.rda")

rm(list=ls())
