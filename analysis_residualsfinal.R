#setwd("~/R/")
#setwd("/Users/jannawilloughby/GDrive/Willoughby lab/furbearer abundance /meta-analysis_ms/Meta analysis/")

library(nlme)
library(ape)
library(phytools)
library(caper)
library(geiger)
library(lme4)
library(scales)

#read in data
data = read.table("alldata_07_29_21.csv", header=T, sep=",")


####1. do non-invasive methods capture the same information as live trapping method?####
#combine camera and genetic into single column - non-invasive methods data
data$noninv = rep(NA, nrow(data))
data$noninv[data$comparison=="livegen"] = data$genetic.num.indv[data$comparison=="livegen"]
data$noninv[data$comparison=="camlive"] = data$camera.num.indv[data$comparison=="camlive"]

#run regression checking if estimates follow 1:1, going through 0,0
lm.out=lm(log(data$noninv)~0+log(data$livetrap.num.indv))
summary(lm.out)
plot(-100,-100, xlab="log(number of individuals), live trapping", ylab="log(number of individuals), non-invasive methods", xlim=c(0,6), ylim=c(0,6))
points(log(data$livetrap.num.indv[data$comparison=="livegen"]), log(data$noninv[data$comparison=="livegen"]), pch=19, col=alpha("firebrick3", alpha=0.7))
points(log(data$livetrap.num.indv[data$comparison=="camlive"]), log(data$noninv[data$comparison=="camlive"]), pch=19, col=alpha("dodgerblue3", alpha=0.7))
segments(x0=0,y0=0, x1=6,y1=(lm.out$coefficients[1]*6), lwd=2, col="black")
data$noninv.resid[data$comparison!="camgen"] = rstandard(lm.out) ### use rstandard() because can interpret compared to 1:1 line, these are standardized residuals from the regression
###Answer: yes, tendency to find more individuals with non invasive compared to live trapping method. How many more?

bootmeans = NULL
alldiffs  = data$noninv-data$livetrap.num.indv
alldiffs  = alldiffs[complete.cases(alldiffs)]
mean(alldiffs)
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 20, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)


####2. do non-invasive methods capture the same information as each other?####
#run regression checking if estimates follow 1:1, going through 0,0
lm.out=lm(log(data$genetic.num.indv)~0+log(data$camera.num.indv))
summary(lm.out)
plot(-100,-100, xlab="log(number of individuals), camera trapping", ylab="log(number of individuals), genetic mark-recapture", xlim=c(0,6), ylim=c(0,6))
points(log(data$camera.num.indv), log(data$genetic.num.indv), pch=19, col=alpha("darkorchid3", alpha=0.7))
segments(x0=0,y0=0, x1=6,y1=(lm.out$coefficients[1]*6), lwd=2, col="black")
data$ncamgen.resid[data$comparison=="camgen"] = rstandard(lm.out)
###Answer: generally more individuals are identified using genetic method compared to camera (points above regression line and line is slope != 1)
bootmeans = NULL
alldiffs  = data$genetic.num.indv-data$camera.num.indv
alldiffs  = alldiffs[complete.cases(alldiffs)]
mean(alldiffs)
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 14, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)

####3. do the differences we observed between trapping/non invasive or within non invasive relate to phylogeny?####
#set up tree
as.phylo.formula2 = function (x, data = parent.frame(), ...){
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if (length(x) != 2) 
    stop(err)
  if (x[[1]] != "~") 
    stop(err)
  f <- x[[2]]
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/") 
      stop(err)
    if (!is.factor(data[[deparse(f[[3]])]])) 
      stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if (length(f) > 1) 
      f <- f[[2]]
  }
  if (!is.factor(data[[deparse(f)]])) 
    stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo)) 
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)])
      t[l] <- paste("(", paste(x, collapse = ","), ")", sep = "")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse = ","),");", sep = "")
  phy <- read.newick(text = string) ## so that singles will be read without error
  phy$edge.length <- rep(1,nrow(phy$edge))
  phy <- collapse.singles(phy)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}
data$order   = as.factor(data$order)
data$family  = as.factor(data$family)
data$genus   = as.factor(data$genus)
data$g_species = as.factor(apply(cbind(as.character(data$genus), as.character(data$species)), 1, paste, collapse="_", sep=""))

#trapping vs non invasive
atree  = as.phylo.formula2(~order/family/genus/g_species, data=data[!is.na(data$noninv.resid),]) 
outpgls = gls(noninv.resid~1, data=data[!is.na(data$noninv.resid),], correlation=corBrownian(1,phy=atree,form=~g_species), method = "ML") #install ape version from pre june 2020? or pgls
summary(outpgls)
###Answer, no. No phylogenetic signal related to residuals.

#compare between camera and genetic
atree  = as.phylo.formula2(~order/family/genus/g_species, data=data[!is.na(data$ncamgen.resid),]) 
outpgls = gls(ncamgen.resid~1, data=data[!is.na(data$ncamgen.resid),], correlation=corBrownian(1,phy=atree,form=~g_species), method = "ML")
summary(outpgls)
###Answer, again no. No phylogenetic signal related to residuals.


####4. do the differences we observed between trapping/non invasive relate technical properties of the different methods?####
##trapping vs noninvasive comparisons
#camera data analysis method
data$camera.method= as.factor(data$camera.method)
lm.out=lm(data$noninv.resid~data$camera.method) 
summary(lm.out)
table(data$camera.method[!is.na(data$noninv.resid) & !is.na(data$camera.method)])
boxplot(data$noninv.resid~data$camera.method, col=c(alpha("dodgerblue3", alpha=0.3), alpha("dodgerblue3", alpha=0.7)))

#bootstrap mean diff in random and spatial comparisons to trapping
tr = data[data$comparison=="camlive" & data$camera.method=="random",]
alldiffs = tr$camera.num.indv-tr$livetrap.num.indv
mean(alldiffs)
bootmeans = NULL
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 3, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)

ts = data[data$comparison=="camlive" & data$camera.method=="spatial",]
alldiffs = ts$camera.num.indv-ts$livetrap.num.indv
mean(alldiffs)
bootmeans = NULL
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 4, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)

#DNA tissue sample collection method
lm.out=lm(data$noninv.resid~data$gen.method)
summary(lm.out)
boxplot(data$noninv.resid~data$gen.method, col=c(alpha("firebrick3", alpha=0.3), alpha("firebrick3", alpha=0.7)))

#bootstrap mean diff in random and spatial comparisons to trapping
tr = data[data$comparison=="camgen" & data$camera.method=="random",]
alldiffs = tr$camera.num.indv-tr$genetic.num.indv
mean(alldiffs)
bootmeans = NULL
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 2, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)

ts = data[data$comparison=="camgen" & data$camera.method=="spatial",]
alldiffs = ts$camera.num.indv-ts$genetic.num.indv
mean(alldiffs)
bootmeans = NULL
for(i in 1:1000){
  bootmeans = c(bootmeans, mean(sample(alldiffs, 12, replace=T)))
}
quantile(bootmeans, probs=c(0.025, 0.975))
mean(bootmeans)


#study site size
lm.out=lm(data$noninv.resid~log(data$study.area.size.km2)) 
summary(lm.out)
plot(-100,-100, xlab="log(study area size) (km2)", ylab="live trapping ~ noninvasive methods, standardized residuals", xlim=c(-1,7), ylim=c(-2,3))
points(log(data$study.area.size.km2), data$noninv.resid, pch=19, col=alpha("darkorchid3", alpha=0.7))
abline(lm.out)

##between noninvasive comparisons
#camera data analysis method
lm.out=lm(data$ncamgen.resid~data$camera.method)
table(data$camera.method[!is.na(data$ncamgen.resid) & !is.na(data$camera.method)])
summary(lm.out)

#DNA tissue sample collection method
lm.out=lm(data$ncamgen.resid~data$gen.method)
summary(lm.out)

#study site size
lm.out=lm(data$ncamgen.resid~log(data$study.area.size.km2)) #this variable has too many levels/most of them in one atm
summary(lm.out)
plot(-100,-100, xlab="log(study area size) (km2)", ylab="genetic mark recapture ~ camera traps, standardized residuals", xlim=c(0,8), ylim=c(-2,3))
points(log(data$study.area.size.km2), data$ncamgen.resid, pch=19, col=alpha("darkorchid3", alpha=0.7))
abline(lm.out)


####5. do the differences we observed relate to body size?####
##trapping vs noninvasive comparison
data$bodysize = apply(cbind(data$bodysize.min.kg, data$bodysize.max.kg), 1, mean, na.rm=T)
lm.out=lm(data$noninv.resid~log(data$bodysize))
summary(lm.out)
plot(-100,-100, xlab="body size (kg)", ylab="live trapping ~ non-invasive methods, standardized residual", xlim=c(-4,5), ylim=c(-2,3))
points(log(data$bodysize[data$comparison=="livegen"]), data$noninv.resid[data$comparison=="livegen"], pch=19, col=alpha("firebrick3", alpha=0.7))
points(log(data$bodysize[data$comparison=="camlive"]), data$noninv.resid[data$comparison=="camlive"], pch=19, col=alpha("dodgerblue3", alpha=0.7))
abline(lm.out)

##compare within noninvasive methods
lm.out=lm(data$ncamgen.resid~log(data$bodysize))
summary(lm.out)
plot(-100,-100, xlab="body size (kg)", ylab="genetic mark-recapture ~ camera traps, standardized residual", xlim=c(-2,6), ylim=c(-1,6))
points(log(data$bodysize[data$comparison=="camgen"]), data$ncamgen.resid[data$comparison=="camgen"], pch=19, col=alpha("darkorchid3", alpha=0.7))
abline(lm.out)



####6. do the differences we observed relate to habitat?####
##trapping vs noninvasive comparison
lm.out=lm(data$noninv.resid~data$habitat.biogeo)
summary(lm.out)
boxplot(data$noninv.resid~data$habita.types.fixed) #differences here, maybe we change this to a half violin plot? boxplots are ugly...
table(data$habita.types.fixed)

##compare within noninvasive methods
lm.out=lm(data$ncamgen.resid~data$habitat.biogeo)
summary(lm.out)
table(data$habita.types.fixed[data$comparison=="camgen"])
boxplot(data$ncamgen.resid~data$habita.types.fixed)



