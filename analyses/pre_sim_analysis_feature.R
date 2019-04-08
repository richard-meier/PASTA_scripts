
library(stats4)
library(rmutil)
options(stringsAsFactors=FALSE)

source("...\\BEINF0_utility_functions.R")

all.met = read.csv("...\\meta_survey_2018_10_10.csv", header=TRUE, check.names=FALSE)

oral.cnt = read.table(
	"...\\rarefied-feature-table-oral.tsv", 
	header=TRUE, skip=1, check.names=FALSE, comment.char = "", sep="\t"
)
colnames(oral.cnt)[1] = "OTU.ID"

tissue.cnt = read.table(
	"...\\rarefied-feature-table-tissue.tsv", 
	header=TRUE, skip=1, check.names=FALSE, comment.char = "", sep="\t"
)
colnames(tissue.cnt)[1] = "OTU.ID"


sbst1 = tissue.cnt[ tissue.cnt$OTU.ID %in% oral.cnt$OTU.ID, ]
sbst2 = tissue.cnt[ !(tissue.cnt$OTU.ID %in% oral.cnt$OTU.ID), ]
cs1 = apply(sbst1[,-1],MARGIN=2,FUN=sum)
cs2 = apply(sbst2[,-1],MARGIN=2,FUN=sum)
sum(cs1)
sum(cs2)
sum(cs1) / (sum(cs1)+sum(cs2))
quantile(cs1/cs2)
plot(cs1,cs2)

sbst1 = oral.cnt[ oral.cnt$OTU.ID %in% tissue.cnt$OTU.ID, ]
sbst2 = oral.cnt[ !(oral.cnt$OTU.ID %in% tissue.cnt$OTU.ID), ]
cs1 = apply(sbst1[,-1],MARGIN=2,FUN=sum)
cs2 = apply(sbst2[,-1],MARGIN=2,FUN=sum)
sum(cs1)
sum(cs2)
sum(cs1) / (sum(cs1)+sum(cs2))
quantile(cs1/cs2)
plot(cs1,cs2)


for(i in 2:length(oral.cnt[1,])){
	oral.cnt[,i] = oral.cnt[,i]/sum(oral.cnt[,i])
}
for(i in 2:length(tissue.cnt[1,])){
	tissue.cnt[,i] = tissue.cnt[,i]/sum(tissue.cnt[,i])
}


exclusionSites = c()
exclusionSites = c( exclusionSites, c("stool","stomach_swab") )

exclusionSet = all.met[all.met$SITE %in% exclusionSites,]$SampleID
all.met = all.met[!(all.met$SITE %in% exclusionSites),]
all.met = all.met[all.met$STUDY != "NDRI",]



consistent_features = intersect( oral.cnt$OTU.ID, tissue.cnt$OTU.ID )
oral.cnt = oral.cnt[oral.cnt$OTU.ID %in% consistent_features,]
tissue.cnt = tissue.cnt[tissue.cnt$OTU.ID %in% consistent_features,]



subo.cnt = oral.cnt
consistentSamples = unique(intersect(colnames(subo.cnt), all.met$SampleID))
subo.cnt = subo.cnt[,c("OTU.ID",consistentSamples)]

subt.cnt = tissue.cnt
consistentSamples = unique(intersect(colnames(subt.cnt), all.met$SampleID))
subt.cnt = subt.cnt[,c("OTU.ID",consistentSamples)]







##############################################################################################################################################################
# TISSUE
sub.cnt = subt.cnt

bla = unlist(sub.cnt[,-1])
bla = bla[bla!=0]
quantile(bla, prob=c(0.05,0.95))
quantile(log(bla), prob=c(0.05,0.95))

OMEGAS_T = c()
PHIS_T = c()
LBOUNDS_T = c()
ALPHAS_T = c()
BETAS_T = c()
SMPLMEANS_T = c()
P0S_T = c()

filepath = "...\\histograms_nonzero_tissue.pdf"

for(i in 1:dim(sub.cnt)[1]){
	ydat = as.matrix(sub.cnt[i,-1])
	P0S_T = c(P0S_T,mean(ydat==0))
	ydat = ydat[ydat!=0]
	NN = length(ydat)
	
	if(NN >= 5){

X = ydat
LL <- function(alpha, beta) {
    return(- sum( dbeta(X,alpha,beta,log=TRUE) ) )
}

mout = mle(LL, start = list(alpha = 1, beta=1))
ALP = coef(mout)[1]
BET = coef(mout)[2]

omega = unname(ALP/(ALP+BET))
phi = unname(ALP + BET)
lbound = 1/max(omega,1-omega)

ALPHAS_T = c(ALPHAS_T,ALP)
BETAS_T = c(BETAS_T,BET)
OMEGAS_T = c(OMEGAS_T,omega)
PHIS_T = c(PHIS_T,phi)
LBOUNDS_T = c(LBOUNDS_T,lbound)
SMPLMEANS_T = c(SMPLMEANS_T,mean(ydat))

	}

}



# process non-zero ######################################################################

bla = unlist(sub.cnt[-1])
bla = bla[bla!=0]
min(bla)

LL <- function(s1, s2) {
    return(- sum( dbeta(OMEGAS_T,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2=1))

plot(density(OMEGAS_T), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="omega")
xx = seq(from=0.001,to=0.2,by=0.001)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(OMEGAS_T,prob=c(0,1/4)),digits=4)
min(OMEGAS_T)/10

omega.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(OMEGAS_T)/10,digits=6)
)



# process 0es ######################################################################

maxprop = 0.999

P0S2 = P0S_T
P0S2[P0S2==1] = maxprop

LL <- function(s1, s2) {
    return(- sum( dbeta(P0S2,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2 = 1))

plot(density(P0S2), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="p0")
xx = seq(from=0.01,to=0.99,by=0.01)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(P0S2,prob=c(0,1/4)),digits=4)
min(P0S2)/10

p0.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(P0S2)/10,digits=6)
)



filepath = "...\\phi_omega_p0_tissue_sub.pdf"
pdf(file=filepath,width=8,height=7)
par(mfrow=c(2,2), mar=c(4.5,5,2,1))
plot(density(OMEGAS_T), xlab=expression(omega) , cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9)
hist(P0S_T,xlab="p : probability of absence", cex.axis=1.2, cex.lab=1.5, main=p0.plotstr, cex.main=0.9)
phifit = lm(log(PHIS_T) ~ log(OMEGAS_T))
omega.plotstr = paste0(
	"b0 = ", round(coef(phifit)[1],4),
	"  and  b1 = ", round(coef(phifit)[2],4)
)
plot(log(OMEGAS_T), log(PHIS_T), xlab=expression(log(omega)), ylab=expression(log(Phi)), cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9  )
xpred = seq(from=-8,to=-1,by=0.1)
lines(xpred,coef(phifit)[1]+xpred*coef(phifit)[2],col="red",lwd=2)
res.str = paste0(
	"mean = ", round(mean(residuals(phifit)),4),
	"  and  sd = ", round(sd(residuals(phifit)),4)
)
plot(density(residuals(phifit)),xlab="residuals", cex=1.5, cex.axis=1.2, cex.lab=1.5, main=res.str, cex.main=0.9 )
dev.off()

round(coef(phifit),digits=2)
round(sd(residuals(phifit)),digits=2)







##############################################################################################################################################################
# ORAL
sub.cnt = subo.cnt

bla = unlist(sub.cnt[,-1])
bla = bla[bla!=0]
quantile(bla, prob=c(0.05,0.95))
quantile(log(bla), prob=c(0.05,0.95))

OMEGAS_O = c()
PHIS_O = c()
LBOUNDS_O = c()
ALPHAS_O = c()
BETAS_O = c()
SMPLMEANS_O = c()
P0S_O = c()

filepath = "...\\histograms_nonzero_tissue.pdf"

for(i in 1:dim(sub.cnt)[1]){
	ydat = as.matrix(sub.cnt[i,-1])
	P0S_O = c(P0S_O,mean(ydat==0))
	ydat = ydat[ydat!=0]
	NN = length(ydat)
	
	if(NN >= 5){

X = ydat
LL <- function(alpha, beta) {
    return(- sum( dbeta(X,alpha,beta,log=TRUE) ) )
}

mout = mle(LL, start = list(alpha = 1, beta=1))
ALP = coef(mout)[1]
BET = coef(mout)[2]

omega = unname(ALP/(ALP+BET))
phi = unname(ALP + BET)
lbound = 1/max(omega,1-omega)

ALPHAS_O = c(ALPHAS_O,ALP)
BETAS_O = c(BETAS_O,BET)
OMEGAS_O = c(OMEGAS_O,omega)
PHIS_O = c(PHIS_O,phi)
LBOUNDS_O = c(LBOUNDS_O,lbound)
SMPLMEANS_O = c(SMPLMEANS_O,mean(ydat))

	}

}



# process non-zero ######################################################################

bla = unlist(sub.cnt[-1])
bla = bla[bla!=0]
min(bla)

LL <- function(s1, s2) {
    return(- sum( dbeta(OMEGAS_O,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2=1))

plot(density(OMEGAS_O), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="omega")
xx = seq(from=0.001,to=0.2,by=0.001)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(OMEGAS_O,prob=c(0,1/4)),digits=4)
min(OMEGAS_O)/10

omega.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(OMEGAS_O)/10,digits=6)
)



# process 0es ######################################################################

maxprop = 0.999

P0S2 = P0S_O
P0S2[P0S2==1] = maxprop

LL <- function(s1, s2) {
    return(- sum( dbeta(P0S2,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2 = 1))

plot(density(P0S2), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="p0")
xx = seq(from=0.01,to=0.99,by=0.01)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(P0S2,prob=c(0,1/4)),digits=4)
min(P0S2)/10

p0.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(P0S2)/10,digits=6)
)



filepath = "...\\phi_omega_p0_oral_sub.pdf"
pdf(file=filepath,width=8,height=7)
par(mfrow=c(2,2), mar=c(4.5,5,2,1))
plot(density(OMEGAS_O), xlab=expression(omega) , cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9)
hist(P0S_O,xlab="p : probability of absence", cex.axis=1.2, cex.lab=1.5, main=p0.plotstr, cex.main=0.9)
phifit = lm(log(PHIS_O) ~ log(OMEGAS_O))
omega.plotstr = paste0(
	"b0 = ", round(coef(phifit)[1],4),
	"  and  b1 = ", round(coef(phifit)[2],4)
)
plot(log(OMEGAS_O), log(PHIS_O), xlab=expression(log(omega)), ylab=expression(log(Phi)), cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9  )
xpred = seq(from=-8,to=-1,by=0.1)
lines(xpred,coef(phifit)[1]+xpred*coef(phifit)[2],col="red",lwd=2)
res.str = paste0(
	"mean = ", round(mean(residuals(phifit)),4),
	"  and  sd = ", round(sd(residuals(phifit)),4)
)
plot(density(residuals(phifit)),xlab="residuals", cex=1.5, cex.axis=1.2, cex.lab=1.5, main=res.str, cex.main=0.9 )
dev.off()







##############################################################################################################################################################
# CONCATENATED

OMEGAS = c(OMEGAS_T,OMEGAS_O)
PHIS = c(PHIS_T,PHIS_O)
LBOUNDS = c(LBOUNDS_T,LBOUNDS_O)
ALPHAS = c(ALPHAS_T,ALPHAS_O)
BETAS = c(BETAS_T,BETAS_O)
SMPLMEANS = c(SMPLMEANS_T,SMPLMEANS_O)
P0S = c(P0S_T,P0S_O)



# process non-zero ######################################################################

LL <- function(s1, s2) {
    return(- sum( dbeta(OMEGAS,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2=1))

plot(density(OMEGAS), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="omega")
xx = seq(from=0.001,to=0.2,by=0.001)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(OMEGAS,prob=c(0,1/4)),digits=4)
min(OMEGAS)/10

omega.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(OMEGAS)/10,digits=6)
)



# process 0es ######################################################################

maxprop = 0.999 #as.numeric(names(tail(table(P0S),n=2))[1])

P0S2 = P0S
P0S2[P0S2==1] = maxprop

LL <- function(s1, s2) {
    return(- sum( dbeta(P0S2,s1,s2,log=TRUE) ) )
}

mout = mle(LL, start = list(s1 = 1, s2 = 1))

plot(density(P0S2), main=paste(c("s1 ="," ; s2 ="),round(coef(mout),digits=2),collapse=""), xlab="p0")
xx = seq(from=0.01,to=0.99,by=0.01)
yy = dbeta(xx,coef(mout)[1],coef(mout)[2])
lines(xx,yy,col="red")

round(coef(mout),digits=2)
round(quantile(P0S2,prob=c(0,1/4)),digits=4)
min(P0S2)/10

p0.plotstr = paste0(
	"(s1, s2) = (",round(coef(mout),digits=2)[1],", ",round(coef(mout),digits=2)[2],")",
	"  and  min/10 = ", round(min(P0S2)/10,digits=6)
)



filepath = "...\\phi_omega_p0_concatenated.pdf"
pdf(file=filepath,width=8,height=7)
par(mfrow=c(2,2), mar=c(4.5,5,2,1))
plot(density(OMEGAS), xlab=expression(omega) , cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9)
hist(P0S,xlab="p : probability of absence", cex.axis=1.2, cex.lab=1.5, main=p0.plotstr, cex.main=0.9)
phifit = lm(log(PHIS) ~ log(OMEGAS))
omega.plotstr = paste0(
	"b0 = ", round(coef(phifit)[1],4),
	"  and  b1 = ", round(coef(phifit)[2],4)
)
plot(log(OMEGAS), log(PHIS), xlab=expression(log(omega)), ylab=expression(log(Phi)), cex=1.5, cex.axis=1.2, cex.lab=1.5, main=omega.plotstr, cex.main=0.9  )
xpred = seq(from=-8,to=-1,by=0.1)
lines(xpred,coef(phifit)[1]+xpred*coef(phifit)[2],col="red",lwd=2)
res.str = paste0(
	"mean = ", round(mean(residuals(phifit)),4),
	"  and  sd = ", round(sd(residuals(phifit)),4)
)
plot(density(residuals(phifit)),xlab="residuals", cex=1.5, cex.axis=1.2, cex.lab=1.5, main=res.str, cex.main=0.9 )
dev.off()

round(coef(phifit),digits=2)
round(sd(residuals(phifit)),digits=2)



filepath = "...\\phi_omega_p0_poster.pdf"
pdf(file=filepath,width=8,height=3.5)
par(mfrow=c(1,2), mar=c(4.5,5,0.2,1.5))
hist(P0S,xlab="p : probability of absence", cex.axis=1.2, cex.lab=1.3, main="", cex.main=0.9,freq=FALSE,breaks=20)
phifit = lm(log(PHIS) ~ log(OMEGAS))
plot(log(OMEGAS), log(PHIS), xlab=expression(log(omega)), ylab=expression(log(Phi)), cex=1.3, cex.axis=1.2, cex.lab=1.3, main="", cex.main=0.9  )
xpred = seq(from=-8,to=-1,by=0.1)
lines(xpred,coef(phifit)[1]+xpred*coef(phifit)[2],col="red",lwd=2)
dev.off()



##############################################################################################################################################################
# compare distributions

x = seq(from = 0.0001, to = 0.1, by = 0.0001)
y1 = dbeta(x=x, shape1=1.46, shape2=121.12)
y2 = dbeta(x=x, shape1=0.82, shape2=59.85)
plot(x,y2,type="l",col="blue", xlab=expression(omega), ylab = expression(f(omega)))
lines(x,y1,col="red")


x = seq(from = 0.001, to = 1, by = 0.001)
y1 = dbeta(x=x, shape1=2.54, shape2=0.53)
y2 = dbeta(x=x, shape1=7.35, shape2=0.49)
plot(x,y2,type="l",col="blue", xlab=expression(omega), ylab = expression(f(omega)))
lines(x,y1,col="red")


x = seq(from = -7, to = -1.5, by = 1)
y1 = -0.89 - 1.1*x
y2 = -1.84 - 1.1*x
plot(x,y2,type="l",col="blue", xlab=expression(omega), ylab = expression(f(omega)))
lines(x,y1,col="red")

