
library(R2OpenBUGS)
MCMC_NUM = 10000
BURN_IN = 1000

options(stringsAsFactors=FALSE)

dataset = read.csv(file="example_data.csv")

nm_y = "relative.abundance"
nm_site = "tissue"
nm_disease = "cancer"
nm_subject = "subject.id"

dataset$GRP_NAME = paste(dataset[,nm_site], dataset[,nm_disease],sep=".")
dataset = dataset[order(dataset$GRP_NAME),]
dataset$GRP = NA
groups = sort(unique(dataset$GRP_NAME))
for(g in 1:length(groups)){
	dataset[dataset$GRP_NAME==groups[g],]$GRP = g
}
dataset$Z = (datatab$Y !=0) * 1.0

atab = dataset[dataset[nm_y]!=0,]
a_subject_grp_tab = data.frame( subject_id=sort(unique(as.vector(atab[,nm_subject]))), a_subject_grp=1:length(unique(atab[,nm_subject])) )
colnames(a_subject_grp_tab)[1] = nm_subject
atab = merge( atab, a_subject_grp_tab, by=nm_subject, all=TRUE )

btab = dataset
btab[nm_y] = (btab[nm_y]==0)*1
b_subject_grp_tab = data.frame( subject_id=sort(unique(as.vector(btab[,nm_subject]))), b_subject_grp=1:length(unique(as.vector(btab[,nm_subject]))) )
colnames(b_subject_grp_tab)[1] = nm_subject
btab = merge( btab, b_subject_grp_tab, by=nm_subject, all=TRUE )

f_mod_logistic <- function() {
    for(i in 1:N){
        Y[i] ~ dbin(p[i], 1)
        logit(p[i]) <- B[GRP[i]] + U[SBJ[i]]
    }
    
    for(j in 1:G){
        B[j] ~ dnorm(0,0.01)
    }
	
	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
	tau ~ dgamma(0.01,0.01)

    phi <- SS*SS
    SS ~ dunif(1,100) 
    
	for(j in 1:gh){
		logit(E_P0_S1[j]) <- B[j]
	}
	
	for(j in ghp1:G){
		logit(E_P0_S2[j]) <- B[j]
	}
    
}

dats = list(
    Y = as.vector(btab[,nm_y]),
    GRP = btab$GRP, 
	SBJ = btab$b_subject_grp, 
	G = max(btab$GRP),
	gh = max(btab$GRP)/2,
	ghp1 = max(btab$GRP)/2 + 1,
	M = length(unique(btab$b_subject_grp)),
    N = length(as.vector((btab[,nm_y])))
)

inivals <- list( list(  B=rep(0,dats$G), tau=2 ) )

fit_logistic <- bugs(
    dats, inits = inivals, parameters.to.save = c("E_P0_S1","E_P0_S2"), 
    model.file = f_mod_logistic, n.chains = 1, n.iter = MCMC_NUM+BURN_IN, n.burnin=BURN_IN, 
    digits=8 #, debug=TRUE
)

P_SITE1 = fit_logistic$sims.list[["E_P0_S1"]]
P_SITE2 = fit_logistic$sims.list[["E_P0_S2"]]

f_mod_beta <- function() {
    for(i in 1:N){
        Y[i] ~ dbeta(a[i], b[i])
        b[i] <- (1-mu[i]) * phi
        a[i] <- mu[i] * phi
        logit(mu[i]) <- B[GRP[i]] + U[SBJ[i]]
    }
	
	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
		
    for(j in 1:G){
        B[j] ~ dnorm(0,0.01)
    }
	tau ~ dgamma(0.01,0.01)
	

	phi <- SS*SS
	SS ~ dunif(1,100)
	
	for(j in 1:gh){
		logit(E_OM_S1[j]) <- B[j]
	}
	
	for(j in ghp1:G){
		logit(E_OM_S2[j]) <- B[j]
	}
	
}

dats = list(
    Y = as.vector(atab[,nm_y]), 
    GRP = atab$GRP, 
	SBJ = atab$a_subject_grp, 
	G = max(atab$GRP),
	gh = max(atab$GRP)/2,
	ghp1 = max(atab$GRP)/2 + 1,
	M = length(unique(atab$a_subject_grp)),
    N = length(as.vector((atab[,nm_y])))
)

inivals <- list( list(  B=rep(0,dats$G), SS=2, tau=2 ) )

fit_beta <- bugs(
    dats, inits=inivals, parameters.to.save = c("B","SS","E_OM_S1","E_OM_S2","phi"), 
    model.file = f_mod_beta, n.chains = 1, n.iter = MCMC_NUM+BURN_IN, n.burnin=BURN_IN, 
    digits=8 #, debug=TRUE
)

O_SITE1 = fit_beta$sims.list[["E_OM_S1"]]
O_SITE2 = fit_beta$sims.list[["E_OM_S2"]]

GH = max(btab$GRP)/2

T3PRSN = c()
T3SPMN = c()
for(j in 1:dim(P_SITE1)[1]){
    T3PRSN[j] = cor( x=P_SITE1[j,1:GH], y=P_SITE2[j,1:GH], method="pearson" )
    T3SPMN[j] = cor( x=P_SITE1[j,1:GH], y=P_SITE2[j,1:GH], method="spearman" )
}

T2PRSN = c()
T2SPMN = c()
for(j in 1:dim(O_SITE1)[1]){
    T2PRSN[j] = cor( x=O_SITE1[j,1:GH], y=O_SITE2[j,1:GH], method="pearson" )
    T2SPMN[j] = cor( x=O_SITE1[j,1:GH], y=O_SITE2[j,1:GH], method="spearman" )
}

muMat = matrix(NA,ncol=max(btab$GRP),nrow=nrow(O_SITE1))
T1PRSN = c()
T1SPMN = c()
for(j in 1:dim(O_SITE1)[1]){
    muMat[j,] = c( O_SITE1[j,1:GH]*(1-P_SITE1[j,1:GH]) , O_SITE2[j,1:GH]*(1-P_SITE2[j,1:GH]) )
    T1PRSN[j] = cor( x=O_SITE1[j,1:GH]*(1-P_SITE1[j,1:GH]), y=O_SITE2[j,1:GH]*(1-P_SITE2[j,1:GH]), method="pearson" )
    T1SPMN[j] = cor( x=O_SITE1[j,1:GH]*(1-P_SITE1[j,1:GH]), y=O_SITE2[j,1:GH]*(1-P_SITE2[j,1:GH]), method="spearman" )
}
omMat = cbind(O_SITE1,O_SITE2)
pMat = cbind(P_SITE1,P_SITE2)

mu_stats = list(
	mean=apply(muMat, MARGIN=2, FUN=mean), 
	lower=apply(muMat, MARGIN=2, FUN=quantile, probs=0.025), 
	upper=apply(muMat, MARGIN=2, FUN=quantile, probs=0.975)
)
omega_stats = list(
	mean=apply(omMat, MARGIN=2, FUN=mean), 
	lower=apply(omMat, MARGIN=2, FUN=quantile, probs=0.025), 
	upper=apply(omMat, MARGIN=2, FUN=quantile, probs=0.975)
)
p_stats = list(
	mean=apply(pMat, MARGIN=2, FUN=mean), 
	lower=apply(pMat, MARGIN=2, FUN=quantile, probs=0.025), 
	upper=apply(pMat, MARGIN=2, FUN=quantile, probs=0.975)
)
model_stats = list(mu=mu_stats,omega=omega_stats,p=p_stats)

TSTATS = cbind(T1PRSN,T2PRSN,T3PRSN,T1SPMN,T2SPMN,T3SPMN)
colnames(TSTATS) = paste0(rep(c("mu.","omega.","p."),2), rep(c("pearson","spearman"),each=3))

h0.prob = function(x, threshold){
	return(mean(x<threshold))
}
p.means = apply(TSTATS, MARGIN=2, FUN=mean)
p.medians = apply(TSTATS, MARGIN=2, FUN=median)
p.lci = apply(TSTATS, MARGIN=2, FUN=quantile,probs=c(0.05))
p.null = apply(TSTATS, MARGIN=2, FUN=h0.prob, threshold=0)

test.results = data.frame(
	parameter = rep(c("mu","omega","p"),2),
	corr.function = rep(c("pearson","spearman"),each=3),
	post.mean = p.means, post.median = p.medians,
	ci.lower = p.lci, h0.probability = p.null
)

par(mfrow=c(3,3))
for(i in 1:ncol(TSTATS)){hist(TSTATS[,i],xlab="T",main=colnames(TSTATS)[i],col="lightgray",xlim=c(-1,1))}

disease_names = sort(unique(dataset[,nm_disease]))
site_names = sort(unique(dataset[,nm_site]))

parnames=c("mu","omega","p")
for(i in 1:length(parnames)){
	nm = parnames[i]
	means = model_stats[[nm]]$mean
	lowers = model_stats[[nm]]$lower
	uppers = model_stats[[nm]]$upper
	plot(
		x=-0.1+(1:GH), y=means[1:GH], ylim=c(0,max(uppers)), xlim=c(1-0.5,GH+0.5), 
		pch=0, col="blue", ylab=nm, xlab="group", main="estimates", xaxt="n"
	)
	axis(side=1, at=(1:GH), labels=disease_names, cex.axis = 1.3)
	for(j in 1:GH){
		lines(x=c(j,j)-0.1,y=c(lowers[j],uppers[j]), col="blue")
	}
	points(x=0.1+(1:GH),y=means[(GH+1):(GH*2)], ylim=c(0,max(uppers)), pch=1, col="red")
	for(j in 1:GH){
		lines(x=c(j,j)+0.1,y=c(lowers[j+GH],uppers[j+GH]), col="red")
	}
}

print(test.results)
