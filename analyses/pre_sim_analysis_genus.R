
library(stats4)
library(rmutil)
options(stringsAsFactors=FALSE)

source("...\\BEINF0_utility_functions.R")

all.met = read.csv("...\\meta_survey_2018_10_10.csv", header=TRUE, check.names=FALSE)
all.cnt = read.csv("...\\final_otu_map_mc100_tax_genus_combined_discrepancy.count.csv", header=TRUE, skip=1, check.names=FALSE)

colnames(all.cnt)[1] = "OTU.ID"

all.met = all.met[!(all.met$subject_id %in% c("378952K","378952")),]

colsums = c()
for(i in 2:length(all.cnt[1,])){
	colsums = c(colsums,sum(all.cnt[,i]))
	all.cnt[,i] = all.cnt[,i]/sum(all.cnt[,i])
}
colsumDat = data.frame(SampleID=colnames(all.cnt)[-1], LIBSIZE=log(colsums))
all.met = merge(all.met,colsumDat, by="SampleID")

# format orst
kingdom = c()
phylum = c()
class = c()
order = c()
family = c()
genus = c()
for( i in 1:length(all.cnt[,1]) ){
	cstr = all.cnt[i,1]
	tax = strsplit(cstr, ";")
	tax = tax[[1]]
	kingdom = c(kingdom, tax[1])
	phylum = c(phylum, tax[2])
	class = c(class, tax[3])
	order = c(order, tax[4])
	family = c(family, tax[5])
	genus = c(genus, tax[6])
}
all.cnt = cbind(kingdom,phylum,class,order,family,genus,all.cnt)

all.cnt[all.cnt$phylum == "p__Saccharibacteria(TM7)",]$phylum = rep("p__Saccharibacteria_(TM7)", sum( all.cnt$phylum == "p__Saccharibacteria(TM7)" ))

exclusionSites = c()
exclusionSites = c( exclusionSites, c("stool","stomach_swab") )

exclusionSet = all.met[all.met$SITE %in% exclusionSites,]$SampleID
all.met = all.met[!(all.met$SITE %in% exclusionSites),]
all.met = all.met[all.met$STUDY != "NDRI",]

all.cnt = all.cnt[,!(colnames(all.cnt) %in% exclusionSet)]
all.met = all.met[!(all.met$SITE %in% exclusionSites),]


sub.cnt = all.cnt
consistentSamples = unique(intersect(colnames(sub.cnt), all.met$SampleID))
sub.cnt = sub.cnt[,c("OTU.ID",consistentSamples)]



bla = unlist(sub.cnt[,-1])
bla = bla[bla!=0]
quantile(bla, prob=c(0.05,0.95))
quantile(log(bla), prob=c(0.05,0.95))

OMEGAS = c()
PHIS = c()
LBOUNDS = c()
ALPHAS = c()
BETAS = c()
SMPLMEANS = c()
P0S = c()

filepath = "...\\histograms_nonzero_tissue.pdf"

for(i in 1:dim(sub.cnt)[1]){
	ydat = as.matrix(sub.cnt[i,-1])
	P0S = c(P0S,mean(ydat==0))
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

ALPHAS = c(ALPHAS,ALP)
BETAS = c(BETAS,BET)
OMEGAS = c(OMEGAS,omega)
PHIS = c(PHIS,phi)
LBOUNDS = c(LBOUNDS,lbound)
SMPLMEANS = c(SMPLMEANS,mean(ydat))

	}

}



# process non-zero ######################################################################

bla = unlist(sub.cnt[-1])
bla = bla[bla!=0]
min(bla)

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

maxprop = 0.999

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



filepath = "...\\phi_omega_p0_genus.pdf"
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
