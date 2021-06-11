## Code for correlating life history with status (threatened/invasive) of legumes
## Updated 08 March 2008 (CJAB)

## Remove everything
rm(list = ls())

## Import data
setwd("~/Documents/Papers/Plants/Legume threat & invasion/data/")
dat <- read.table("legume4.csv",header=T,sep=",")

## Import NatureServe data for NA species only
#dat <- read.table("C:\\Documents and Settings\\cbradshaw\\My Documents\\Papers\\Extinctions\\Legume threat & invasion\\data\\natureservelegumes.csv",header=T,sep=",")


## Test response
#test.var <- "threat"
test.var <- "IUCN"
#test.var <- "invasive"
#test.var <- "both"

## univariate models only
univ <- "no"
#univ <- "yes"

## Remove IUCN threat classes B and D2 (i.e., those 'restricted in range')
rng.rem <- 1 ## '1' to remove B&D2 species, '0' to keep them

## With or without 'seedsize'
#wseedsize <- 1
wseedsize <- 0

## With height as covariate
#wheightcov <- 1
wheightcov <- 0

## With leaf length as covariate
#wllmaxcov <- 1
wllmaxcov <- 0
 
## With nominate phylogenetic groups (NPG) instead of subfamily?
wNPG <- "no" ## "no" if random effect should be subfam/tribe

 
#Convert var to factors
dat$subfam <- factor(dat$subfam)
dat$tribe <- factor(dat$tribe)
#dat$ecores <- factor(ecores)
dat$ndist <- factor(dat$ndist)
dat$curdist <- factor(dat$curdist)
dat$NPG <- factor(dat$NPG)
#dat$range <- factor(ifelse(range == 2 | range == 3, 2, 1))
dat$range <- factor(dat$range, ordered=T)
dat$endem <- factor(dat$endem)
dat$alt <- factor(dat$alt)
dat$hbt <- factor(dat$hbt)
dat$cult <- factor(dat$cult)
#dat$hab <- factor(dat$hab)
dat$hab <- factor(dat$hab, levels = rev(sort(unique(dat$hab))),ordered=T)
dat$mode <- factor(dat$mode)
dat$ht <- factor(dat$ht, ordered=T)
lhtcov <- log10(ifelse(dat$ht == 1, 0.5, (ifelse(dat$ht == 2, 5, 15))))
dat$arm <- factor(dat$arm)
dat$hair <- factor(dat$hair)
dat$leaflmin <- factor(dat$leaflmin)
#dat$leaflmax <- factor(dat$leaflmax)
dat$leaflmax <- factor(dat$leaflmax,ordered=T)
lllmaxcov <- log10(dat$llmaxcov)
dat$leafmarg <- factor(dat$leafmarg)
dat$leaftype <- factor(dat$leaftype)
dat$flordisp <- factor(dat$flordisp)
dat$hook <- factor(dat$hook)
dat$deh <- factor(dat$deh)
dat$seedsize <- factor(dat$seedsize)
dat$seedcoat <- factor(dat$seedcoat)
dat <- data.frame(dat,lhtcov,lllmaxcov)

## Examine NA numbers (%)
attach(dat)
100*(1-(length(na.omit(as.numeric(range)))/length(as.numeric(range))))
100*(1-(length(na.omit(as.numeric(endem)))/length(as.numeric(endem))))
100*(1-(length(na.omit((alt)))/length((alt))))
100*(1-(length(na.omit((hbt)))/length((hbt))))
100*(1-(length(na.omit((cult)))/length((cult))))
100*(1-(length(na.omit((hab)))/length((hab))))
100*(1-(length(na.omit((mode)))/length((mode))))
100*(1-(length(na.omit((ht)))/length((ht))))
100*(1-(length(na.omit((arm)))/length((arm))))
100*(1-(length(na.omit((hair)))/length((hair))))
100*(1-(length(na.omit((leaflmax)))/length((leaflmax))))
100*(1-(length(na.omit((leafmarg)))/length((leafmarg))))
100*(1-(length(na.omit((leaftype)))/length((leaftype))))
100*(1-(length(na.omit((flordisp)))/length((flordisp))))
100*(1-(length(na.omit((hook)))/length((hook))))
100*(1-(length(na.omit((deh)))/length((deh))))
100*(1-(length(na.omit((seedsize)))/length((seedsize))))
100*(1-(length(na.omit((seedcoat)))/length((seedcoat))))
detach(dat)

## Correlation matrix
dat.cor <- na.omit(data.frame(as.integer(dat$range),as.integer(dat$alt),as.integer(dat$hbt),as.integer(dat$hab),as.integer(dat$mode),as.integer(dat$ht),as.integer(dat$arm),as.integer(dat$hair),as.integer(dat$leaflmax),as.integer(dat$flordisp),as.integer(dat$hook),as.integer(dat$deh),as.integer(dat$seedsize)))
colnames(dat.cor) <- c("range","altitude","habitat","habit","life cycle","height","armaments","hair","leaf length","fdisplay","hooks","deh","ssize")
cor(dat.cor,y=NULL,method="spearman")

cor(dat.cor,y=NULL,method="kendall")


## IUCN categories
#IUCN status 
#1=EX or EW; 2=CR; 3=EN; 4=VU; 5=LR-cd/nt (NT); 6=DD (removed); 7=Normal; 8=Weedy; 9=Invasive 
summary(as.factor(dat$IUCN))
IUCNstat <- as.factor(ifelse(dat$IUCN <= 4, 1, 0))
dat <- data.frame(dat,IUCNstat)

## Choose data subset
## Remove invasive species for threat/not threat analysis
if (test.var == "threat") {
data <- dat[dat$ecores != 1,]

## Change threat code
threat <- as.factor(abs(as.integer(data$ecores)))
data <- data.frame(data,threat)}

## Remove invasive, weedy and normal species for IUCN comparison (EX,EW,CR,EN,VU versus LR)
if (test.var == "IUCN") {
data <- dat[dat$IUCN <=5,]}

## Remove threatened species for invasive/not invasive analysis
if (test.var == "invasive") {
data <- dat[dat$ecores != -1,]
inv <- as.factor(abs(as.integer(data$ecores)))
data <- data.frame(data,inv)}

## Remove 'normal' species for invasive/threatened analysis
if (test.var == "both") {
data <- dat[dat$ecores != 0,]
thinv <- as.factor(ifelse(as.integer(data$ecores) == -1, 1, 0))
data <- data.frame(data,thinv)}

## Remove range-restricted species (IUCN class B2 & D) (for threat)
iucnrge <- ifelse(data$IUCNtyp == "B" | data$IUCNtyp == "D2" | data$IUCNtyp == "BD2", 0, 1)
data <- data.frame(data,iucnrge)
if (rng.rem == 1) {
data <- data[data$iucnrge == 1,]}

## remove NAs for threat status analysis
if (test.var == "threat" & wseedsize == 1) {
data1 <- na.omit(data.frame(data$cult,data$threat,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$alt,data$hbt,data$hab,data$mode,data$ht,data$arm,data$hair,data$leaflmax,data$flordisp,data$hook,data$deh,data$seedsize))
colnames(data1) <- c("cult","threat","subfam","NPG","tribe","gen","range","alt","hbt","hab","mode","ht","arm","hair","leaflmax","flordisp","hook","deh","seedsize")}

if (test.var == "threat" & wseedsize == 0) {
data1 <- na.omit(data.frame(data$cult,data$threat,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$alt,data$hbt,data$hab,data$mode,data$ht,data$arm,data$hair,data$leaflmax,data$flordisp,data$hook,data$deh))
colnames(data1) <- c("cult","threat","subfam","NPG","tribe","gen","range","alt","hbt","hab","mode","ht","arm","hair","leaflmax","flordisp","hook","deh")}

if (test.var == "threat" & wseedsize == 1 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$cult,data$threat,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$alt,data$hbt,data$hab,data$mode,data$lhtcov,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$hook,data$deh,data$seedsize))
colnames(data1) <- c("cult","threat","subfam","NPG","tribe","gen","range","alt","hbt","hab","mode","ht","arm","hair","leaflmax","flordisp","hook","deh","seedsize")}

if (test.var == "threat" & wseedsize == 0 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$cult,data$threat,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$alt,data$hbt,data$hab,data$mode,data$lhtcov,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$hook,data$deh))
colnames(data1) <- c("cult","threat","subfam","NPG","tribe","gen","range","alt","hbt","hab","mode","ht","arm","hair","leaflmax","flordisp","hook","deh")}


## remove NAs for IUCN status analysis
if (test.var == "IUCN") {
data1 <- na.omit(data.frame(data$cult,data$IUCNstat,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$alt,data$hbt,data$hab,data$mode,data$ht,data$arm,data$hair,data$leaflmax,data$flordisp,data$hook,data$deh))
colnames(data1) <- c("cult","IUCNstat","subfam","NPG","tribe","gen","range","alt","hbt","hab","mode","ht","arm","hair","leaflmax","flordisp","hook","deh")}


## remove NAs for invasiveness analysis
if (test.var == "invasive" & wseedsize == 1) {
data1 <- na.omit(data.frame(data$inv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$leaflmax,data$flordisp,data$ht,data$deh,data$hook,data$alt,data$hbt,data$seedsize))
colnames(data1) <- c("inv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt","seedsize")}

if (test.var == "invasive" & wseedsize == 0) {
data1 <- na.omit(data.frame(data$inv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$leaflmax,data$flordisp,data$ht,data$deh,data$hook,data$alt,data$hbt))
colnames(data1) <- c("inv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt")}

if (test.var == "invasive" & wseedsize == 1 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$inv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$lhtcov,data$deh,data$hook,data$alt,data$hbt,data$seedsize))
colnames(data1) <- c("inv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt","seedsize")}

if (test.var == "invasive" & wseedsize == 0 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$inv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$lhtcov,data$deh,data$hook,data$alt,data$hbt))
colnames(data1) <- c("inv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt")}


## remove NAs for both analysis
if (test.var == "both" & wseedsize == 1) {
data1 <- na.omit(data.frame(data$thinv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$leaflmax,data$flordisp,data$ht,data$deh,data$hook,data$alt,data$hbt,data$seedsize))
colnames(data1) <- c("thinv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt","seedsize")}

if (test.var == "both" & wseedsize == 0) {
data1 <- na.omit(data.frame(data$thinv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$leaflmax,data$flordisp,data$ht,data$deh,data$hook,data$alt,data$hbt))
colnames(data1) <- c("thinv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt")}

if (test.var == "both" & wseedsize == 1 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$thinv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$lhtcov,data$deh,data$hook,data$alt,data$hbt,data$seedsize))
colnames(data1) <- c("thinv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt","seedsize")}

if (test.var == "both" & wseedsize == 0 & wheightcov == 1 & wllmaxcov == 1) {
data1 <- na.omit(data.frame(data$thinv,data$subfam,data$NPG,data$tribe,data$gen,data$range,data$hab,data$mode,data$arm,data$hair,data$lllmaxcov,data$flordisp,data$lhtcov,data$deh,data$hook,data$alt,data$hbt))
colnames(data1) <- c("thinv","subfam","NPG","tribe","gen","range","hab","mode","arm","hair","leaflmax","flordisp","ht","deh","hook","alt","hbt")}


## remove cultivars (for threat & IUCN)
if (test.var == "threat" | test.var == "IUCN") {
data1 <- data1[data1$cult == 0,]}

## Sample size
dim(data1)[1]
if (test.var == "invasive") {
  table(data1$inv)}
  
if (test.var == "threat") {
  table(data1$threat)}

if (test.var == "IUCN") {
  table(data1$IUCNstat)}

if (test.var == "both") {
  table(data1$thinv)}

## select model type
model.type <- "lmer"
#model.type <- "glm"

## Load libraries
library(MASS)
library(lme4)
library(Matrix)

AICc.lmer <- function(...) {
 models <- list(...)
 num.mod <- length(models)
 AICcs <- numeric(num.mod)
 ns <- numeric(num.mod)
 ks <- numeric(num.mod)
 AICc.vec <- rep(0,num.mod)
 for (i in 1:num.mod) {
  n <- models[[i]]@XtX@x[1]
  k <- length(models[[i]]@fixef)+sum(models[[i]]@nc)+1
  AICcs[i] <- (-2*logLik(models[[i]])[1]) + ((2*k*n)/(n-k-1))
  ns[i] <- n
  ks[i] <- k
  AICc.vec[i] <- AICcs[i]}
 return(AICc.vec)}

BIC.lmer <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	BICs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	BIC.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
    n <- models[[i]]@XtX@x[1]
    k <- length(models[[i]]@fixef)+sum(models[[i]]@nc)+1
		BICs[i] <- (-2*logLik(models[[i]])) + k*log(n)
		ns[i] <- n
		ks[i] <- k
		BIC.vec[i] <- BICs[i]}
	return(BIC.vec)}

# Set functions
AICc.glm <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	AICcs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	AICc.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
		if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
		if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
		AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
		ns[i] <- n
		ks[i] <- k
		AICc.vec[i] <- AICcs[i]}
	return(AICc.vec)}

BIC.glm <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	BICs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	BIC.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
		if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
		if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
		BICs[i] <- (-2*logLik(models[[i]])) + k*log(n)
		ns[i] <- n
		ks[i] <- k
		BIC.vec[i] <- BICs[i]}
	return(BIC.vec)}
                 
k.lmer <- function(x) {
  if (length(x$df.residual) == 0) length(x@fixef)+sum(x@nc)+1}

k.glm <- function(x) {
  if (length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)}
 
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
chdev.lmer <- function(x) ((( as.numeric(dev.null) - as.numeric(deviance(x)))/ as.numeric(dev.null)*100))
chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object

n <- dim(data1)[1]
n

#######################
## Define model sets ##
#######################

##########
## lmer ##
##########

## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "threat" & wseedsize == 1 & univ == "no") {

## Height-only model
mod.1 <- "threat~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "threat~range+(1|subfam/tribe)"
mod.3 <- "threat~range+hbt+(1|subfam/tribe)"
mod.4 <- "threat~range+alt+(1|subfam/tribe)"
mod.5 <- "threat~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "threat~mode+(1|subfam/tribe)"
mod.7 <- "threat~hab+(1|subfam/tribe)"
mod.8 <- "threat~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "threat~flordisp+(1|subfam/tribe)"
mod.10 <- "threat~hook+(1|subfam/tribe)"
mod.11 <- "threat~deh+(1|subfam/tribe)"
mod.12 <- "threat~seedsize+(1|subfam/tribe)"
mod.13 <- "threat~hook+deh+(1|subfam/tribe)"
mod.14 <- "threat~flordisp+hook+deh+(1|subfam/tribe)"
mod.15 <- "threat~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "threat~arm+(1|subfam/tribe)"
mod.17 <- "threat~hair+(1|subfam/tribe)"
mod.18 <- "threat~leaflmax+(1|subfam/tribe)"
mod.19 <- "threat~arm+hair+(1|subfam/tribe)"
mod.20 <- "threat~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "threat~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "threat~range+hab+(1|subfam/tribe)"
mod.23 <- "threat~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "threat~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "threat~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "threat~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "threat~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "threat~deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.29 <- "threat~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "threat~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "threat~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.32 <- "threat~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "threat~ht+range+(1|subfam/tribe)"
mod.34 <- "threat~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "threat~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "threat~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "threat~ht+mode+(1|subfam/tribe)"
mod.38 <- "threat~ht+hab+(1|subfam/tribe)"
mod.39 <- "threat~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "threat~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "threat~ht+hook+(1|subfam/tribe)"
mod.42 <- "threat~ht+deh+(1|subfam/tribe)"
mod.43 <- "threat~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "threat~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "threat~ht+flordisp+hook+deh+(1|subfam/tribe)"
mod.46 <- "threat~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "threat~ht+arm+(1|subfam/tribe)"
mod.48 <- "threat~ht+hair+(1|subfam/tribe)"
mod.49 <- "threat~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "threat~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "threat~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "threat~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "threat~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "threat~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "threat~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "threat~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "threat~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "threat~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "threat~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.60 <- "threat~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "threat~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "threat~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.63 <- "threat~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "threat~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "threat~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13,mod.14,mod.15,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.29,mod.30,mod.31,mod.32,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.43,mod.44,mod.45,mod.46,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.60,mod.61,mod.62,mod.63,mod.64,mod.65)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "threat" & univ == "yes") {

mod.1 <- "threat~ht+(1|subfam/tribe)"
mod.2 <- "threat~range+(1|subfam/tribe)"
mod.3 <- "threat~hbt+(1|subfam/tribe)"
mod.4 <- "threat~alt+(1|subfam/tribe)"
mod.5 <- "threat~mode+(1|subfam/tribe)"
mod.6 <- "threat~hab+(1|subfam/tribe)"
mod.7 <- "threat~flordisp+(1|subfam/tribe)"
mod.8 <- "threat~hook+(1|subfam/tribe)"
mod.9 <- "threat~deh+(1|subfam/tribe)"
mod.10 <- "threat~arm+(1|subfam/tribe)"
mod.11 <- "threat~hair+(1|subfam/tribe)"
mod.12 <- "threat~leaflmax+(1|subfam/tribe)"
mod.13 <- "threat~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13)
}




## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "threat" & wseedsize == 0 & univ == "no") {

## Height-only model
mod.1 <- "threat~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "threat~range+(1|subfam/tribe)"
mod.3 <- "threat~range+hbt+(1|subfam/tribe)"
mod.4 <- "threat~range+alt+(1|subfam/tribe)"
mod.5 <- "threat~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "threat~mode+(1|subfam/tribe)"
mod.7 <- "threat~hab+(1|subfam/tribe)"
mod.8 <- "threat~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "threat~flordisp+(1|subfam/tribe)"
mod.10 <- "threat~hook+(1|subfam/tribe)"
mod.11 <- "threat~deh+(1|subfam/tribe)"
#mod.12 <- "threat~seedsize+(1|subfam/tribe)"
mod.13 <- "threat~hook+deh+(1|subfam/tribe)"
mod.14 <- "threat~flordisp+hook+deh+(1|subfam/tribe)"
#mod.15 <- "threat~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "threat~arm+(1|subfam/tribe)"
mod.17 <- "threat~hair+(1|subfam/tribe)"
mod.18 <- "threat~leaflmax+(1|subfam/tribe)"
mod.19 <- "threat~arm+hair+(1|subfam/tribe)"
mod.20 <- "threat~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "threat~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "threat~range+hab+(1|subfam/tribe)"
mod.23 <- "threat~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "threat~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "threat~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "threat~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "threat~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "threat~deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.29 <- "threat~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "threat~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "threat~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.32 <- "threat~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "threat~ht+range+(1|subfam/tribe)"
mod.34 <- "threat~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "threat~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "threat~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "threat~ht+mode+(1|subfam/tribe)"
mod.38 <- "threat~ht+hab+(1|subfam/tribe)"
mod.39 <- "threat~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "threat~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "threat~ht+hook+(1|subfam/tribe)"
mod.42 <- "threat~ht+deh+(1|subfam/tribe)"
#mod.43 <- "threat~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "threat~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "threat~ht+flordisp+hook+deh+(1|subfam/tribe)"
#mod.46 <- "threat~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "threat~ht+arm+(1|subfam/tribe)"
mod.48 <- "threat~ht+hair+(1|subfam/tribe)"
mod.49 <- "threat~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "threat~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "threat~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "threat~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "threat~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "threat~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "threat~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "threat~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "threat~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "threat~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "threat~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.60 <- "threat~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "threat~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "threat~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.63 <- "threat~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "threat~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "threat~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.13,mod.14,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.30,mod.31,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.44,mod.45,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.61,mod.62,mod.64,mod.65)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "threat" & wseedsize == 0 & univ == "no" & wNPG == "yes") {

## Height-only model
mod.1 <- "threat~ht+(1|NPG/tribe)"

## 'habitat' models
mod.2 <- "threat~range+(1|NPG/tribe)"
mod.3 <- "threat~range+hbt+(1|NPG/tribe)"
mod.4 <- "threat~range+alt+(1|NPG/tribe)"
mod.5 <- "threat~range+hbt+alt+(1|NPG/tribe)"

## 'life history' models
mod.6 <- "threat~mode+(1|NPG/tribe)"
mod.7 <- "threat~hab+(1|NPG/tribe)"
mod.8 <- "threat~mode+hab+(1|NPG/tribe)"

## 'reproduction' models
mod.9 <- "threat~flordisp+(1|NPG/tribe)"
mod.10 <- "threat~hook+(1|NPG/tribe)"
mod.11 <- "threat~deh+(1|NPG/tribe)"
#mod.12 <- "threat~seedsize+(1|NPG/tribe)"
mod.13 <- "threat~hook+deh+(1|NPG/tribe)"
mod.14 <- "threat~flordisp+hook+deh+(1|NPG/tribe)"
#mod.15 <- "threat~flordisp+hook+deh+seedsize+(1|NPG/tribe)"

## 'defence' models
mod.16 <- "threat~arm+(1|NPG/tribe)"
mod.17 <- "threat~hair+(1|NPG/tribe)"
mod.18 <- "threat~leaflmax+(1|NPG/tribe)"
mod.19 <- "threat~arm+hair+(1|NPG/tribe)"
mod.20 <- "threat~leaflmax+arm+hair+(1|NPG/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "threat~range+mode+hab+(1|NPG/tribe)"
mod.22 <- "threat~range+hab+(1|NPG/tribe)"
mod.23 <- "threat~range+hbt+mode+hab+(1|NPG/tribe)"
mod.24 <- "threat~range+alt+mode+hab+(1|NPG/tribe)"
mod.25 <- "threat~range+hbt+alt+mode+hab+(1|NPG/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "threat~flordisp+arm+hair+leaflmax+(1|NPG/tribe)"
mod.27 <- "threat~hook+arm+hair+leaflmax+(1|NPG/tribe)"
mod.28 <- "threat~deh+arm+hair+leaflmax+(1|NPG/tribe)"
#mod.29 <- "threat~seedsize+arm+hair+leaflmax+(1|NPG/tribe)"
mod.30 <- "threat~hook+deh+arm+hair+leaflmax+(1|NPG/tribe)"
mod.31 <- "threat~flordisp+hook+deh+arm+hair+leaflmax+(1|NPG/tribe)"
#mod.32 <- "threat~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|NPG/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "threat~ht+range+(1|NPG/tribe)"
mod.34 <- "threat~ht+range+hbt+(1|NPG/tribe)"
mod.35 <- "threat~ht+range+alt+(1|NPG/tribe)"
mod.36 <- "threat~ht+range+hbt+alt+(1|NPG/tribe)"

## 'life history' models
mod.37 <- "threat~ht+mode+(1|NPG/tribe)"
mod.38 <- "threat~ht+hab+(1|NPG/tribe)"
mod.39 <- "threat~ht+mode+hab+(1|NPG/tribe)"

## 'reproduction' models
mod.40 <- "threat~ht+flordisp+(1|NPG/tribe)"
mod.41 <- "threat~ht+hook+(1|NPG/tribe)"
mod.42 <- "threat~ht+deh+(1|NPG/tribe)"
#mod.43 <- "threat~ht+seedsize+(1|NPG/tribe)"
mod.44 <- "threat~ht+hook+deh+(1|NPG/tribe)"
mod.45 <- "threat~ht+flordisp+hook+deh+(1|NPG/tribe)"
#mod.46 <- "threat~ht+flordisp+hook+deh+seedsize+(1|NPG/tribe)"

## 'defence' models
mod.47 <- "threat~ht+arm+(1|NPG/tribe)"
mod.48 <- "threat~ht+hair+(1|NPG/tribe)"
mod.49 <- "threat~ht+leaflmax+(1|NPG/tribe)"
mod.50 <- "threat~ht+arm+hair+(1|NPG/tribe)"
mod.51 <- "threat~ht+leaflmax+arm+hair+(1|NPG/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "threat~ht+range+mode+hab+(1|NPG/tribe)"
mod.53 <- "threat~ht+range+hab+(1|NPG/tribe)"
mod.54 <- "threat~ht+range+hbt+mode+hab+(1|NPG/tribe)"
mod.55 <- "threat~ht+range+alt+mode+hab+(1|NPG/tribe)"
mod.56 <- "threat~ht+range+hbt+alt+mode+hab+(1|NPG/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "threat~ht+flordisp+arm+hair+leaflmax+(1|NPG/tribe)"
mod.58 <- "threat~ht+hook+arm+hair+leaflmax+(1|NPG/tribe)"
mod.59 <- "threat~ht+deh+arm+hair+leaflmax+(1|NPG/tribe)"
#mod.60 <- "threat~ht+seedsize+arm+hair+leaflmax+(1|NPG/tribe)"
mod.61 <- "threat~ht+hook+deh+arm+hair+leaflmax+(1|NPG/tribe)"
mod.62 <- "threat~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|NPG/tribe)"
#mod.63 <- "threat~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|NPG/tribe)"

## saturated model
mod.64 <- "threat~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax+(1|NPG/tribe)"

## null model
mod.65 <- "threat~1+(1|NPG/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.13,mod.14,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.30,mod.31,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.44,mod.45,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.61,mod.62,mod.64,mod.65)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "IUCN" & wseedsize == 0 & univ == "no") {

## Height-only model
mod.1 <- "IUCNstat~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "IUCNstat~range+(1|subfam/tribe)"
mod.3 <- "IUCNstat~hbt+(1|subfam/tribe)"
mod.4 <- "IUCNstat~alt+(1|subfam/tribe)"

## 'life history' models
mod.5 <- "IUCNstat~mode+(1|subfam/tribe)"
mod.6 <- "IUCNstat~hab+(1|subfam/tribe)"

## 'reproduction' models
mod.7 <- "IUCNstat~flordisp+(1|subfam/tribe)"
mod.8 <- "IUCNstat~hook+(1|subfam/tribe)"
mod.9 <- "IUCNstat~deh+(1|subfam/tribe)"

## 'defence' models
mod.10 <- "IUCNstat~arm+(1|subfam/tribe)"
mod.11<- "IUCNstat~hair+(1|subfam/tribe)"
mod.12 <- "IUCNstat~leaflmax+(1|subfam/tribe)"

## null model
mod.13 <- "IUCNstat~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.13)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "invasive" & wseedsize == 1 & univ == "no") {

## Height-only model
mod.1 <- "inv~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "inv~range+(1|subfam/tribe)"
mod.3 <- "inv~range+hbt+(1|subfam/tribe)"
mod.4 <- "inv~range+alt+(1|subfam/tribe)"
mod.5 <- "inv~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "inv~mode+(1|subfam/tribe)"
mod.7 <- "inv~hab+(1|subfam/tribe)"
mod.8 <- "inv~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "inv~flordisp+(1|subfam/tribe)"
mod.10 <- "inv~hook+(1|subfam/tribe)"
mod.11 <- "inv~deh+(1|subfam/tribe)"
mod.12 <- "inv~seedsize+(1|subfam/tribe)"
mod.13 <- "inv~hook+deh+(1|subfam/tribe)"
mod.14 <- "inv~flordisp+hook+deh+(1|subfam/tribe)"
mod.15 <- "inv~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "inv~arm+(1|subfam/tribe)"
mod.17 <- "inv~hair+(1|subfam/tribe)"
mod.18 <- "inv~leaflmax+(1|subfam/tribe)"
mod.19 <- "inv~arm+hair+(1|subfam/tribe)"
mod.20 <- "inv~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "inv~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "inv~range+hab+(1|subfam/tribe)"
mod.23 <- "inv~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "inv~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "inv~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "inv~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "inv~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "inv~deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.29 <- "inv~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "inv~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "inv~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.32 <- "inv~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "inv~ht+range+(1|subfam/tribe)"
mod.34 <- "inv~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "inv~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "inv~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "inv~ht+mode+(1|subfam/tribe)"
mod.38 <- "inv~ht+hab+(1|subfam/tribe)"
mod.39 <- "inv~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "inv~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "inv~ht+hook+(1|subfam/tribe)"
mod.42 <- "inv~ht+deh+(1|subfam/tribe)"
mod.43 <- "inv~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "inv~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "inv~ht+flordisp+hook+deh+(1|subfam/tribe)"
mod.46 <- "inv~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "inv~ht+arm+(1|subfam/tribe)"
mod.48 <- "inv~ht+hair+(1|subfam/tribe)"
mod.49 <- "inv~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "inv~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "inv~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "inv~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "inv~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "inv~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "inv~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "inv~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "inv~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "inv~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "inv~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.60 <- "inv~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "inv~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "inv~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.63 <- "inv~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "inv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "inv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13,mod.14,mod.15,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.29,mod.30,mod.31,mod.32,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.43,mod.44,mod.45,mod.46,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.60,mod.61,mod.62,mod.63,mod.64,mod.65)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "invasive" & univ == "yes") {

mod.1 <- "inv~ht+(1|subfam/tribe)"
mod.2 <- "inv~range+(1|subfam/tribe)"
mod.3 <- "inv~hbt+(1|subfam/tribe)"
mod.4 <- "inv~alt+(1|subfam/tribe)"
mod.5 <- "inv~mode+(1|subfam/tribe)"
mod.6 <- "inv~hab+(1|subfam/tribe)"
mod.7 <- "inv~flordisp+(1|subfam/tribe)"
mod.8 <- "inv~hook+(1|subfam/tribe)"
mod.9 <- "inv~deh+(1|subfam/tribe)"
mod.10 <- "inv~arm+(1|subfam/tribe)"
mod.11 <- "inv~hair+(1|subfam/tribe)"
mod.12 <- "inv~leaflmax+(1|subfam/tribe)"
mod.13 <- "inv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13)
}



## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "invasive" & wseedsize == 0 & univ == "no") {

## Height-only model
mod.1 <- "inv~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "inv~range+(1|subfam/tribe)"
mod.3 <- "inv~range+hbt+(1|subfam/tribe)"
mod.4 <- "inv~range+alt+(1|subfam/tribe)"
mod.5 <- "inv~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "inv~mode+(1|subfam/tribe)"
mod.7 <- "inv~hab+(1|subfam/tribe)"
mod.8 <- "inv~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "inv~flordisp+(1|subfam/tribe)"
mod.10 <- "inv~hook+(1|subfam/tribe)"
mod.11 <- "inv~deh+(1|subfam/tribe)"
#mod.12 <- "inv~seedsize+(1|subfam/tribe)"
mod.13 <- "inv~hook+deh+(1|subfam/tribe)"
mod.14 <- "inv~flordisp+hook+deh+(1|subfam/tribe)"
#mod.15 <- "inv~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "inv~arm+(1|subfam/tribe)"
mod.17 <- "inv~hair+(1|subfam/tribe)"
mod.18 <- "inv~leaflmax+(1|subfam/tribe)"
mod.19 <- "inv~arm+hair+(1|subfam/tribe)"
mod.20 <- "inv~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "inv~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "inv~range+hab+(1|subfam/tribe)"
mod.23 <- "inv~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "inv~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "inv~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "inv~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "inv~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "inv~deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.29 <- "inv~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "inv~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "inv~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.32 <- "inv~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "inv~ht+range+(1|subfam/tribe)"
mod.34 <- "inv~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "inv~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "inv~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "inv~ht+mode+(1|subfam/tribe)"
mod.38 <- "inv~ht+hab+(1|subfam/tribe)"
mod.39 <- "inv~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "inv~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "inv~ht+hook+(1|subfam/tribe)"
mod.42 <- "inv~ht+deh+(1|subfam/tribe)"
#mod.43 <- "inv~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "inv~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "inv~ht+flordisp+hook+deh+(1|subfam/tribe)"
#mod.46 <- "inv~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "inv~ht+arm+(1|subfam/tribe)"
mod.48 <- "inv~ht+hair+(1|subfam/tribe)"
mod.49 <- "inv~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "inv~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "inv~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "inv~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "inv~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "inv~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "inv~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "inv~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "inv~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "inv~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "inv~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.60 <- "inv~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "inv~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "inv~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.63 <- "inv~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "inv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "inv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.13,mod.14,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.30,mod.31,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.44,mod.45,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.61,mod.62,mod.64,mod.65)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "both" & wseedsize == 1 & univ == "no") {

## Height-only model
mod.1 <- "thinv~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "thinv~range+(1|subfam/tribe)"
mod.3 <- "thinv~range+hbt+(1|subfam/tribe)"
mod.4 <- "thinv~range+alt+(1|subfam/tribe)"
mod.5 <- "thinv~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "thinv~mode+(1|subfam/tribe)"
mod.7 <- "thinv~hab+(1|subfam/tribe)"
mod.8 <- "thinv~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "thinv~flordisp+(1|subfam/tribe)"
mod.10 <- "thinv~hook+(1|subfam/tribe)"
mod.11 <- "thinv~deh+(1|subfam/tribe)"
mod.12 <- "thinv~seedsize+(1|subfam/tribe)"
mod.13 <- "thinv~hook+deh+(1|subfam/tribe)"
mod.14 <- "thinv~flordisp+hook+deh+(1|subfam/tribe)"
mod.15 <- "thinv~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "thinv~arm+(1|subfam/tribe)"
mod.17 <- "thinv~hair+(1|subfam/tribe)"
mod.18 <- "thinv~leaflmax+(1|subfam/tribe)"
mod.19 <- "thinv~arm+hair+(1|subfam/tribe)"
mod.20 <- "thinv~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "thinv~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "thinv~range+hab+(1|subfam/tribe)"
mod.23 <- "thinv~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "thinv~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "thinv~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "thinv~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "thinv~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "thinv~deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.29 <- "thinv~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "thinv~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "thinv~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.32 <- "thinv~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "thinv~ht+range+(1|subfam/tribe)"
mod.34 <- "thinv~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "thinv~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "thinv~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "thinv~ht+mode+(1|subfam/tribe)"
mod.38 <- "thinv~ht+hab+(1|subfam/tribe)"
mod.39 <- "thinv~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "thinv~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "thinv~ht+hook+(1|subfam/tribe)"
mod.42 <- "thinv~ht+deh+(1|subfam/tribe)"
mod.43 <- "thinv~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "thinv~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "thinv~ht+flordisp+hook+deh+(1|subfam/tribe)"
mod.46 <- "thinv~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "thinv~ht+arm+(1|subfam/tribe)"
mod.48 <- "thinv~ht+hair+(1|subfam/tribe)"
mod.49 <- "thinv~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "thinv~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "thinv~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "thinv~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "thinv~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "thinv~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "thinv~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "thinv~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "thinv~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "thinv~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "thinv~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.60 <- "thinv~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "thinv~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "thinv~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.63 <- "thinv~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "thinv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "thinv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13,mod.14,mod.15,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.29,mod.30,mod.31,mod.32,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.43,mod.44,mod.45,mod.46,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.60,mod.61,mod.62,mod.63,mod.64,mod.65)
}



## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "both" & univ == "yes") {

mod.1 <- "thinv~ht+(1|subfam/tribe)"
mod.2 <- "thinv~range+(1|subfam/tribe)"
mod.3 <- "thinv~hbt+(1|subfam/tribe)"
mod.4 <- "thinv~alt+(1|subfam/tribe)"
mod.5 <- "thinv~mode+(1|subfam/tribe)"
mod.6 <- "thinv~hab+(1|subfam/tribe)"
mod.7 <- "thinv~flordisp+(1|subfam/tribe)"
mod.8 <- "thinv~hook+(1|subfam/tribe)"
mod.9 <- "thinv~deh+(1|subfam/tribe)"
mod.10 <- "thinv~arm+(1|subfam/tribe)"
mod.11 <- "thinv~hair+(1|subfam/tribe)"
mod.12 <- "thinv~leaflmax+(1|subfam/tribe)"
mod.13 <- "thinv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.12,mod.13)
}


## New a priori mod.vec sets
if (model.type == "lmer" & test.var == "both" & wseedsize == 0 & univ == "no") {

## Height-only model
mod.1 <- "thinv~ht+(1|subfam/tribe)"

## 'habitat' models
mod.2 <- "thinv~range+(1|subfam/tribe)"
mod.3 <- "thinv~range+hbt+(1|subfam/tribe)"
mod.4 <- "thinv~range+alt+(1|subfam/tribe)"
mod.5 <- "thinv~range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.6 <- "thinv~mode+(1|subfam/tribe)"
mod.7 <- "thinv~hab+(1|subfam/tribe)"
mod.8 <- "thinv~mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.9 <- "thinv~flordisp+(1|subfam/tribe)"
mod.10 <- "thinv~hook+(1|subfam/tribe)"
mod.11 <- "thinv~deh+(1|subfam/tribe)"
#mod.12 <- "thinv~seedsize+(1|subfam/tribe)"
mod.13 <- "thinv~hook+deh+(1|subfam/tribe)"
mod.14 <- "thinv~flordisp+hook+deh+(1|subfam/tribe)"
#mod.15 <- "thinv~flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.16 <- "thinv~arm+(1|subfam/tribe)"
mod.17 <- "thinv~hair+(1|subfam/tribe)"
mod.18 <- "thinv~leaflmax+(1|subfam/tribe)"
mod.19 <- "thinv~arm+hair+(1|subfam/tribe)"
mod.20 <- "thinv~leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.21 <- "thinv~range+mode+hab+(1|subfam/tribe)"
mod.22 <- "thinv~range+hab+(1|subfam/tribe)"
mod.23 <- "thinv~range+hbt+mode+hab+(1|subfam/tribe)"
mod.24 <- "thinv~range+alt+mode+hab+(1|subfam/tribe)"
mod.25 <- "thinv~range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.26 <- "thinv~flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.27 <- "thinv~hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.28 <- "thinv~deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.29 <- "thinv~seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.30 <- "thinv~hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.31 <- "thinv~flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.32 <- "thinv~flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## All above with allometric control (height)
## 'habitat' models
mod.33 <- "thinv~ht+range+(1|subfam/tribe)"
mod.34 <- "thinv~ht+range+hbt+(1|subfam/tribe)"
mod.35 <- "thinv~ht+range+alt+(1|subfam/tribe)"
mod.36 <- "thinv~ht+range+hbt+alt+(1|subfam/tribe)"

## 'life history' models
mod.37 <- "thinv~ht+mode+(1|subfam/tribe)"
mod.38 <- "thinv~ht+hab+(1|subfam/tribe)"
mod.39 <- "thinv~ht+mode+hab+(1|subfam/tribe)"

## 'reproduction' models
mod.40 <- "thinv~ht+flordisp+(1|subfam/tribe)"
mod.41 <- "thinv~ht+hook+(1|subfam/tribe)"
mod.42 <- "thinv~ht+deh+(1|subfam/tribe)"
#mod.43 <- "thinv~ht+seedsize+(1|subfam/tribe)"
mod.44 <- "thinv~ht+hook+deh+(1|subfam/tribe)"
mod.45 <- "thinv~ht+flordisp+hook+deh+(1|subfam/tribe)"
#mod.46 <- "thinv~ht+flordisp+hook+deh+seedsize+(1|subfam/tribe)"

## 'defence' models
mod.47 <- "thinv~ht+arm+(1|subfam/tribe)"
mod.48 <- "thinv~ht+hair+(1|subfam/tribe)"
mod.49 <- "thinv~ht+leaflmax+(1|subfam/tribe)"
mod.50 <- "thinv~ht+arm+hair+(1|subfam/tribe)"
mod.51 <- "thinv~ht+leaflmax+arm+hair+(1|subfam/tribe)"

## 'habitat' + 'life history' models combined
mod.52 <- "thinv~ht+range+mode+hab+(1|subfam/tribe)"
mod.53 <- "thinv~ht+range+hab+(1|subfam/tribe)"
mod.54 <- "thinv~ht+range+hbt+mode+hab+(1|subfam/tribe)"
mod.55 <- "thinv~ht+range+alt+mode+hab+(1|subfam/tribe)"
mod.56 <- "thinv~ht+range+hbt+alt+mode+hab+(1|subfam/tribe)"

## 'reproduction' + 'defence' & 'leaf type' models combined
mod.57 <- "thinv~ht+flordisp+arm+hair+leaflmax+(1|subfam/tribe)"
mod.58 <- "thinv~ht+hook+arm+hair+leaflmax+(1|subfam/tribe)"
mod.59 <- "thinv~ht+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.60 <- "thinv~ht+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"
mod.61 <- "thinv~ht+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
mod.62 <- "thinv~ht+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"
#mod.63 <- "thinv~ht+flordisp+hook+deh+seedsize+arm+hair+leaflmax+(1|subfam/tribe)"

## saturated model
mod.64 <- "thinv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax+(1|subfam/tribe)"

## null model
mod.65 <- "thinv~1+(1|subfam/tribe)"

## Make model vector
mod.vec <- c(mod.1,mod.2,mod.3,mod.4,mod.5,mod.6,mod.7,mod.8,mod.9,mod.10,mod.11,mod.13,mod.14,mod.16,mod.17,mod.18,mod.19,mod.20,mod.21,mod.22,mod.23,mod.24,mod.25,mod.26,mod.27,mod.28,mod.30,mod.31,mod.33,mod.34,mod.35,mod.36,mod.37,mod.38,mod.39,mod.40,mod.41,mod.42,mod.44,mod.45,mod.47,mod.48,mod.49,mod.50,mod.51,mod.52,mod.53,mod.54,mod.55,mod.56,mod.57,mod.58,mod.59,mod.61,mod.62,mod.64,mod.65)
}

##########
## glm ##
##########


## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)

  if (model.type == "lmer" & test.var == "threat") dev.null <- deviance(lmer(threat~1+(1|subfam/tribe),family=binomial(link="logit"),data=data1,method="Laplace",control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000))))
  if (model.type == "lmer" & test.var == "threat" & wNPG == "yes") dev.null <- deviance(lmer(threat~1+(1|NPG/tribe),family=binomial(link="logit"),data=data1,method="Laplace",control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000))))
  if (model.type == "lmer" & test.var == "IUCN") dev.null <- deviance(lmer(IUCNstat~1+(1|subfam/tribe),family=binomial(link="logit"),data=data1,method="Laplace",control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000))))
  if (model.type == "lmer" & test.var == "invasive") dev.null <- deviance(lmer(inv~1+(1|subfam/tribe),family=binomial(link="logit"),data=data1,method="Laplace",control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000))))
  if (model.type == "lmer" & test.var == "both") dev.null <- deviance(lmer(thinv~1+(1|subfam/tribe),family=binomial(link="logit"),data=data1,method="Laplace",control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000))))

  mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {

    if (model.type == "lmer") fit <- lmer(as.formula(mod.vec[i]),family=binomial(link="logit"), data=data1, method="Laplace", na.action=na.omit,control=(list(PQLmaxIt=1000,maxIter=1000,msMaxIter=1000,niterEM=1000)))
  	if (model.type == "glm") fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), data=data1)

    LL.vec[i] <- as.numeric(logLik(fit))
    if (model.type == "lmer") AICc.vec[i] <- AICc.lmer(fit)
    if (model.type == "glm") AICc.vec[i] <- AICc.glm(fit)
    BIC.vec[i] <- ifelse(model.type == "lmer",BIC.lmer(fit),BIC.glm(fit))
    k.vec[i] <- ifelse(model.type == "lmer",k.lmer(fit),k.glm(fit))
  	if (model.type == "lmer") pc.dev.vec[i] <- chdev.lmer(fit)
    if (model.type == "glm") pc.dev.vec[i] <- chdev.glm(fit)

  	print(i)
}

#AIC weights
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)

#AIC weights
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

## Create results dataframe
table<-cbind(mod.num,k.vec,LL.vec,AICc.vec,dAICc,wAICc,BIC.vec,dBIC,wBIC,pc.dev.vec)
colnames(table)<-c("model","k","-LogL","AICc","dAICc","wAIC","BIC","dBIC","wBIC","pcdev")
rownames(table)<- mod.vec

##shows table sorted by wAIC
summary.table<-table[order(table[,6],decreasing=TRUE),1:10]
summary.table[1:10,]

##shows table sorted by wBIC
summary.table<-table[order(table[,9],decreasing=TRUE),1:10]
summary.table[1:10,]


## Sample size
dim(data1)[1]
if (test.var == "invasive") {
  table(data1$inv)}
  
if (test.var == "threat") {
  table(data1$threat)}

if (test.var == "IUCN") {
  table(data1$IUCNstat)}

if (test.var == "both") {
  table(data1$thinv)}



##################
## Predictions
##################

## Test response
#test.var <- "threat"
#test.var <- "invasive"
test.var <- "both"

## Based on saturated model?
#mod.sat <- "sat"
mod.sat <- "no"

## need to use glmmPQL to use 'predict' function
if (test.var == "threat") fit.pred <- glmmPQL(threat~ht+range+hbt+alt+mode+hab,random = ~1|subfam/tribe, family=binomial(link=logit), data=data1)
if (test.var == "invasive") fit.pred <- glmmPQL(inv~ht+range+hbt+alt+mode+hab,random = ~1|subfam/tribe, family=binomial(link=logit), data=data1)
if (test.var == "both") fit.pred <- glmmPQL(thinv~ht+range+hbt+alt+mode+hab,random = ~1|subfam/tribe, family=binomial(link=logit), data=data1)

if (test.var == "invasive" & mod.sat == "sat") fit.pred <- glmmPQL(inv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax,random = ~1|subfam/tribe, family=binomial(link=logit), data=data1)
if (test.var == "both" & mod.sat == "sat") fit.pred <- glmmPQL(thinv~ht+range+hbt+alt+mode+hab+flordisp+hook+deh+arm+hair+leaflmax,random = ~1|subfam/tribe, family=binomial(link=logit), data=data1)


rm(data.test)
data.test <- data1
## Height data (1 = <=1 m, 2 = 1 to <=10 m, 3 = >10 m)
#ht.new <- factor(rep(3,dim(data.test)[1]))
#data.test$ht <- ht.new

## Range data (1 = 1 kingdom, 2 = 2 kingdoms, 3 = 3 kingdoms)
#range.new <- factor(rep(3,dim(data.test)[1]))
#data.test$range <- range.new

## Habitat data (1 = closed forest, 2 = open habitat; 3 = both)
#hbt.new <- factor(rep(3,dim(data.test)[1]))
#data.test$hbt <- hbt.new

## Altitude data (1 = lowland, 2 = montane; 3 = both)
#alt.new <- factor(rep(3,dim(data.test)[1]))
#data.test$alt <- alt.new

## Life cycle data (1 = annual, 2 = perennial)
#mode.new <- factor(rep(2,dim(data.test)[1]))
#data.test$mode <- mode.new

## Habit data (1 = tree, 2 = shrub, 3 = climber, 4 = herb)
hab.new <- factor(rep(4,dim(data.test)[1]))
data.test$hab <- hab.new

## Floral display data (1 = solitary, 2 = inflorescence)
#flordisp.new <- factor(rep(2,dim(data.test)[1]))
#data.test$flordisp <- flordisp.new

## Fruit hook data (1 = present, 0 = absent)
#hook.new <- factor(rep(0,dim(data.test)[1]))
#data.test$hook <- hook.new

## Fruit dehiscence data (1 = present, 0 = absent)
#deh.new <- factor(rep(1,dim(data.test)[1]))
#data.test$deh <- deh.new

## Armaments data (1 = present, 0 = absent)
#arm.new <- factor(rep(0,dim(data.test)[1]))
#data.test$arm <- arm.new

## Vegetative hair data (1 = present, 0 = absent)
#hair.new <- factor(rep(0,dim(data.test)[1]))
#data.test$hair <- hair.new

## Max leaf lamina length data (1 = < 4 cm; 2 = 4-8 cm, 3 = 8-12 cm, 4 = 12-16 cm, 5 = 16-20 cm, 6 = >20 cm)
#leaflmax.new <- factor(rep(6,dim(data.test)[1]))
#data.test$leaflmax <- leaflmax.new


pred.out <- na.omit(predict(fit.pred,data.test,type="response"))
pred.base <- na.omit(predict(fit.pred,data1,type="response"))
iter <- 10000
mean.out <- base.out <- rep(0,iter)
	for (i in 1:iter) {
		resamp1 <- sample(pred.out,replace=T)
    resamp0 <- sample(pred.base,replace=T)
		mean.out[i] <- mean(resamp1); base.out[i] <- mean(resamp0)
		#print(i)
	} ## end i loop
mean.up95 <- (quantile(mean.out,probs=0.975))
mean.up95
mean.lo95 <- (quantile(mean.out,probs=0.025))
mean(na.omit(predict(fit.pred,data.test,type="response")))
mean.lo95

## bootstrap predicted from original data
meanbase.up95 <- (quantile(base.out,probs=0.975))
meanbase.up95
meanbase.lo95 <- (quantile(base.out,probs=0.025))
mean(na.omit(predict(fit.pred,data1,type="response")))
meanbase.lo95
