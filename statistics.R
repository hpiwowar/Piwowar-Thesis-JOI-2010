#install.packages('Hmisc')
#install.packages('Design')
#install.packages('tree')
#install.packages('tree')

library(Hmisc,T)
library(Design,T)
library(tree)

#### READ DATA
setwd(".")
setwd("Documents/Projects/ASIST Metrics 2009/JofInformetrics/")
dat.raw = read.csv("rawdata.csv", header=TRUE, sep=",")
dim(dat.raw)
names(dat.raw) = gsub("_", ".", names(dat.raw))
names(dat.raw)


### SET UP DERIVED VARIABLES
first.career.length = 2008-dat.raw["first.first.year"]
names(first.career.length) = c("first.career.length")
last.career.length = 2008-dat.raw["last.first.year"]
names(last.career.length) = c("last.career.length")

dat.orig = cbind(dat.raw, first.career.length, last.career.length)
dat.orig$policy.strength = ordered(dat.orig$policy.strength)


#### Compute principal components of author experience
### FIRST AUTHOR
pc.first = princomp(scale( cbind(log(1+dat.orig$first.hindex), 
								log(1+dat.orig$first.aindex), 
								dat.orig$first.career.length)))
summary(pc.first)
# NOTE THE NEGATIVE!  This is so that
# higher hindexes etc correlate with higher author experience scores
first.author.exp = - pc.first$scores[,1]  
pc.first$loadings
#biplot(pc.first)
first.author.exp.freq = ecdf(first.author.exp)(first.author.exp)

### LAST AUTHOR
pc.last = princomp(scale( cbind(log(1+dat.orig$last.hindex), 
								log(1+dat.orig$last.aindex), 
								dat.orig$last.career.length)))
summary(pc.last)
# NOTE THE NEGATIVE!  This is so that
# higher hindexes etc correlate with higher author experience scores
last.author.exp = - pc.last$scores[,1]
pc.last$loadings
#biplot(pc.last)
last.author.exp.freq = ecdf(last.author.exp)(last.author.exp)

#### OVERVIEW TABLE

v.orig = names(dat.orig)[3:length(dat.orig)]
v = c("is.data.shared", "policy.strength", "impact.factor", "is.usa.address", 
"is.nih.funded", "any.nih.data.sharing.applies", "sum.of.max.award.for.each.grant", 
"any.direct.cost.over.500k", "any.new.or.renewed.since.2003", "num.authors",
"first.author.exp", "last.author.exp", 
"first.career.length", "last.career.length"
)
dat.extra = cbind(dat.orig, first.author.exp, last.author.exp)
dat = dat.extra[,v]						
dd = datadist(dat)
options(datadist='dd')
options(digits=2)

describe(dat.orig[,v.orig], listunique=0) 
summary(dat.orig[,v.orig]) 
	
##### FIGURE 1  (mean line was added manually later)

cuteq = function(X, n) cut(X,quantile(X,(0:n)/n),include.lowest=TRUE) 
	
dots = NULL
dots$is.data.shared = dat.orig$is.data.shared
dots$impact.factor = cut(dat.orig$impact.factor, c(0,6,8,15,50)) 
dots$journal.policy.strength = dat.orig$policy.strength
dots$num.authors = cuteq(dat.orig$num.authors,3)
dots$FIRST.author.hindex = cuteq(dat.orig$first.hindex,3) 
dots$FIRST.author.aindex = cuteq(dat.orig$first.aindex,3) 
dots$FIRST.author.career.length = cuteq(dat.orig$first.career.length,3) 
dots$LAST.author.hindex = cuteq(dat.orig$last.hindex,3) 
dots$LAST.author.aindex = cuteq(dat.orig$last.aindex,3)
dots$LAST.author.career.length = cuteq(dat.orig$last.career.length,3) 
dots$is.usa.address = dat.orig$is.usa.address
dots$is.nih.funded = dat.orig$is.nih.funded
dots$nih.funds = cut(dat.orig$sum.of.max.award.for.each.grant/1000, c(-1, 1, 750, 2000, 500000))
dots$nih.requires.data.sharing.plan = dat.orig$any.nih.data.sharing.applies
dots$is.nih.grant.number.missing = is.na(dots$nih.requires.data.sharing.plan)

s = summary(is.data.shared ~ ., dat=dots)
s


tiff("figure1_no_mean_line.tiff", bg="white", width=880, height=1200)
plot(s)
title("Proportion of studies with shared datasets")
dev.off()


### DO IMPUTATION

do.imputation = function(column) {
	dat.orig$temp = na.include(column)
	mytree = tree(temp ~ policy.strength + num.authors + 
		is.usa.address + 
		is.nih.funded + 
		first.hindex + first.aindex + first.career.length + 
		first.num.papers + first.total.pmc.citations + 
		last.hindex + last.aindex + last.career.length + 
		last.num.papers + last.total.pmc.citations +
		log(impact.factor),
		dat=dat.orig, 
		control=tree.control(nobs=502, mincut=20) 
		)
	summary(mytree)
	plot(mytree); text(mytree)
	response = impute(column, 
			predict(mytree, dat.orig)[is.na(column)])
	response = round(response)
	print(response)
	response
}

par(mfrow=c(1,2))
dat$any.nih.data.sharing.applies.imputed = 
	do.imputation(dat.orig$any.nih.data.sharing.applies)
dat$sum.of.max.award.for.each.grant.imputed = 
	do.imputation(dat.orig$sum.of.max.award.for.each.grant)

#######  MULTIVARIATE REGRESSION

par(mfrow=c(1,1))
dat.all = cbind(abs(dots$is.nih.grant.number.missing), dat.extra, dat[15:18])

dat$log.award = log(1+sum.of.max.award.for.each.grant.imputed)

dd = datadist(dat)
options(datadist='dd')
options(digits=2)
	
f = lrm(formula = is.data.shared ~ ordered(policy.strength) + 
	is.usa.address*is.nih.funded + 
	any.nih.data.sharing.applies.imputed + 
	rcs(first.author.exp,4) + rcs(last.author.exp, 4) +
	rcs(log(1+sum.of.max.award.for.each.grant.imputed), 4) + 
	rcs(log(impact.factor), 4),
	dat=dat.all, x=T, y=T
	)
anova(f)
f
par(mfrow=c(4,5))
resid(f, 'partial', pl=TRUE)
resid(f, 'gof')

#### TABLE 1
anova(f)


#####  MODEL FOR ODDS RATIO ANALYSIS

f2 = lrm(formula = is.data.shared ~ ordered(policy.strength) + 
	is.usa.address + 
	any.nih.data.sharing.applies.imputed + 
	rcs(first.author.exp,4) + rcs(last.author.exp, 4) +
	rcs(log(1+sum.of.max.award.for.each.grant.imputed), 4) + 
	rcs(log(impact.factor), 4),
	dat=dat, x=T, y=T
	)
f2
anova(f2)
par(mfrow=c(4,5))
resid(f2, 'partial', pl=TRUE)
resid(f2, 'gof')

# The p-value is large (>.3 in that example) indicating no significant lack of fit.
# from http://www.unc.edu/courses/2006spring/ecol/145/001/docs/solutions/final.htm
residuals.lrm(f2,type='gof')


####  FIGURE 2
attach(dat)

tiff("figure2.tiff", bg="white", width=880)

summ = summary(f2, 
	sum.of.max.award.for.each.grant.imputed=c(1,750000),
	policy.strength=0,
	impact.factor=c(5,15), 
	first.author.exp=quantile(first.author.exp, c(0.25,0.75)),
	last.author.exp=quantile(last.author.exp, c(0.25,0.75))
	)
summ
plot(summ, log=T)

dev.off()

### FIGURE 3

tiff("figure3.tiff", bg="white", width=880)

cutHML = function(X) {
	n=3
    cuts = cut(X,quantile(X,(0:n)/n),include.lowest=TRUE, labels=FALSE)
	ifelse(cuts==1, 'low', ifelse(cuts==2, 'med', 'high'))
	}

par(mfrow=c(1,2))  
group = factor(cutHML(first.author.exp))
grouplty=NULL
grouplty[levels(group)] = c(1:3)
labels = c("high", "med", "low")
plsmo(log(impact.factor), is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	xlim=c(1,4), ylim=c(0,1), label.curves=F,
	ylab="Probability of shared data", xlab= "Log of impact factor")
legend("topleft", labels, lty=grouplty[labels], bty="n")
title("IF by First author experience")

group = factor(cutHML(last.author.exp))
grouplty=NULL
grouplty[levels(group)] = c(1:3)
labels = c("high", "med", "low")
plsmo(log(impact.factor), is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	xlim=c(1,4), ylim=c(0,1), label.curves=F,
	ylab = "Probability of shared data", 
	xlab="Log of impact factor")
legend("topleft", labels, lty=grouplty[labels], bty="n")	
title("IF by Last author experience")	

dev.off()


#### FIGURE 4

tiff("figure4.tiff", bg="white", width=880)

par(mfrow=c(1,3))  
group = policy.strength
levels(group) = c("None", "Weak", "Strong")  #### Hrm, not very robust!
grouplabels = c("Strong", "Weak", "None")
grouplty=NULL
grouplty[grouplabels] = c(3:1)


plsmo(log(impact.factor), is.data.shared, group=group, lty=grouplty,
	datadensity=F,
	xlim=c(1,4), ylim=c(0,1), trim=0, label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "Log of impact factor")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")	
title("IF by Journal policy strength")			
plsmo(first.author.exp, is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	ylim=c(0,1), label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "First author experience")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")		
title("First author by Journal policy")				
plsmo(last.author.exp, is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	ylim=c(0,1), label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "Last author experience")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")			
title("Last author by Journal policy")

dev.off()

##### FIGURE 5

tiff("figure5.tiff", bg="white", width=880)

par(mfrow=c(1,3))  
group = factor(any.nih.data.sharing.applies)
#### Hrm, this isn't a very robust ordering of the legend.
## It is accurate for this version of the data, but may not be accurate if the data are resorted.
## What is a better way to do this?
levels(group) = c("Not required", "NIH sharing plan required")  
grouplabels = c("NIH sharing plan required", "Not required")
grouplty=NULL
grouplty[grouplabels] = c(2:1)

plsmo(log(impact.factor), is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	xlim=c(1,4), ylim=c(0,1), trim=0, label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "Log of impact factor")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")			
title("IF by NIH policy")	
plsmo(first.author.exp, is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	ylim=c(0,1), label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "First author experience")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")			
title("First author by NIH policy")	
plsmo(last.author.exp, is.data.shared, group=group, lty=grouplty,
	datadensity=F, 
	ylim=c(0,1), label.curves=F,
	ylab = "Probability of shared data", 
	xlab= "Last author experience")
legend("topleft", grouplabels, lty=grouplty[levels(group)], bty="n")			
title("Last author by NIH policy")	

dev.off()
