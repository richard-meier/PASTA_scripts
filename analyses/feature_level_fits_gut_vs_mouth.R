
library(R2OpenBUGS)
library(mcmcse)
library(coda)

bugs.script.hacked <-
  function(parameters.to.save, n.chains, n.iter, n.burnin,
           n.thin, saveExec, restart, model.file.bug,
           model.file, debug=FALSE, is.inits, 
           DIC=FALSE, useWINE=FALSE,
           newWINE=TRUE, WINEPATH=NULL, bugs.seed=NULL, summary.only=FALSE,
           save.history=(.Platform$OS.type == "windows" | useWINE==TRUE),
           bugs.data.file, bugs.inits.files,
           over.relax = FALSE)
{
    print("Utilizing hacked openbugs script to increase chain precision!")
  ## Write file script.txt for Bugs
  if(n.iter - n.burnin < 2)
    stop ("(n.iter-n.burnin) must be at least 2")
  working.directory <- getwd()
  script <- "script.txt"

  model <- 
    if (length(grep("\\\\", model.file)) || length(grep("/", model.file))) {
        gsub("\\\\", "/", model.file)
    }
    else file.path(working.directory, model.file)
  model <- R2OpenBUGS:::native2win(model, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  data <- file.path(working.directory, bugs.data.file)
  data <- R2OpenBUGS:::native2win(data, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  coda  <- file.path(working.directory, "/")
  coda <- R2OpenBUGS:::native2win(coda, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  model.file.bug<-file.path(working.directory,model.file.bug)
  model.file.bug<-R2OpenBUGS:::native2win(model.file.bug, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  logFile <- file.path(working.directory, "log.odc")
  logFile <- R2OpenBUGS:::native2win(logFile, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)
  logFileTxt <- file.path(working.directory, "log.txt")
  logFileTxt <- R2OpenBUGS:::native2win(logFileTxt, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  inits <- paste(working.directory, "/", bugs.inits.files, sep="")
  inits <- sapply(inits, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH, 
    function(x, useWINE, newWINE, WINEPATH) 
    {R2OpenBUGS:::native2win(x, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)})

  initlist <- paste("modelInits(", "'", inits, "',",1:n.chains,")\n", sep="")

  savelist <- paste("samplesSet(", parameters.to.save, ")\n", sep="")
  summarylist <- paste("summarySet(", parameters.to.save, ")\n", sep="")

  bugs.seed.cmd <- ""
  if (!is.null(bugs.seed)) {
        bugs.seed.cmd <- paste("modelSetRN(", bugs.seed, ")\n", sep="")
  }
    
  thinUpdate <- paste("modelUpdate(", formatC(n.burnin, format='d'), ",", n.thin, 
                      ",",formatC(n.burnin, format='d'), ")\n", sep="")

  cat(
    if(.Platform$OS.type == "windows" | useWINE) "modelDisplay('log')\n",
    if(restart)c("modelInternalize('", model.file.bug, "')\n"),
    if(restart && n.burnin>0)c(
       "samplesClear('*')\n",
       "summaryClear('*')\n"
       ),
    "modelPrecision(6)\n",
    if(!restart)c( 
      "modelCheck('", model, "')\n",
      "modelData('", data, "')\n",
      "modelCompile(", n.chains, ")\n"
      ),
    if(!restart)bugs.seed.cmd,
    if(!restart && is.inits) initlist,
    if(!restart)"modelGenInits()\n",
    if(!restart && over.relax) 'over.relax("yes")\n',
    if((!restart) || (n.burnin>0))c(
    thinUpdate,
    savelist,
    summarylist
    ),
    if(((!restart) || (n.burnin>0)) && DIC) "dicSet()\n",
    "modelUpdate(", formatC(n.iter-n.burnin, format='d'), ",", n.thin, 
                  ",",formatC(n.iter-n.burnin, format='d'),")\n", 
    "samplesCoda('*', '", coda, "')\n", 
    "summaryStats('*')\n", 
    if(DIC) "dicStats()\n",
    if (save.history) "samplesHistory('*')\n", 
    if(saveExec)c("modelExternalize('",model.file.bug,"')\n"),
    if(.Platform$OS.type == "windows" | useWINE) c("modelSaveLog('", logFile, "')\n",
    "modelSaveLog('", logFileTxt, "')\n"),
    file=script, sep="", append=FALSE)

    if(!debug) cat("modelQuit('y')\n", file=script, append=TRUE)

  sims.files <- paste("CODAchain", 1:n.chains, ".txt", sep="")
  for(i in 1:n.chains)
    cat("OpenBUGS did not run correctly.\n", file=sims.files[i], append=FALSE)
}

assignInNamespace("bugs.script", bugs.script.hacked, ns="R2OpenBUGS")










########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
#######
#######    gut vs mouth general + random effects (3 ICDC)
#######
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

options(stringsAsFactors=FALSE)
source("...\\BEINF0_utility_functions.R")
MCMC_NUM = 5000
ALPHA = 0.1

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

colsums_o = c()
for(i in 2:length(oral.cnt[1,])){
	colsums_o = c(colsums_o,sum(oral.cnt[,i]))
	oral.cnt[,i] = oral.cnt[,i]/sum(oral.cnt[,i])
}
colsumDat_o = data.frame(SampleID=colnames(oral.cnt)[-1], LIBSIZE=log(colsums_o))

colsums_t = c()
for(i in 2:length(tissue.cnt[1,])){
	colsums_t = c(colsums_t,sum(tissue.cnt[,i]))
	tissue.cnt[,i] = tissue.cnt[,i]/sum(tissue.cnt[,i])
}
colsumDat_t = data.frame(SampleID=colnames(tissue.cnt)[-1], LIBSIZE=log(colsums_t))

colsumDat = rbind(colsumDat_o,colsumDat_t)
all.met = merge(all.met,colsumDat,all=TRUE,by=c("SampleID"))



exclusionSites = c()
exclusionSites = c( exclusionSites, c("stool","stomach_swab","intest_duodenal_stent", "panc_swab") )

exclusionSet = all.met[all.met$SITE %in% exclusionSites,]$SampleID
all.met = all.met[!(all.met$SITE %in% exclusionSites),]
all.met = all.met[all.met$STUDY != "NDRI",]
all.met = all.met[!(all.met$SampleID %in% c("491762.1","491762.3","491762.4","491762.5")),]



consistent_features = intersect( oral.cnt$OTU.ID, tissue.cnt$OTU.ID )
oral.cnt = oral.cnt[oral.cnt$OTU.ID %in% consistent_features,]
tissue.cnt = tissue.cnt[tissue.cnt$OTU.ID %in% consistent_features,]



consistentSamples = unique(
	union(
		intersect(colnames(oral.cnt), all.met$SampleID),
		intersect(colnames(tissue.cnt), all.met$SampleID)
	)
)
all.met = all.met[all.met$SampleID %in% consistentSamples,]
subo.cnt = oral.cnt
subo.cnt = subo.cnt[,c("OTU.ID",colnames(subo.cnt)[colnames(subo.cnt) %in% consistentSamples])]
subt.cnt = tissue.cnt
subt.cnt = subt.cnt[,c("OTU.ID",colnames(subt.cnt)[colnames(subt.cnt) %in% consistentSamples])]

sub.cnt = merge(subt.cnt,subo.cnt,all=TRUE,by="OTU.ID")



MEDAGE = median(all.met$AGE)
MEDBMI = median(na.omit(all.met$BMI))

newGender = c()
for(i in 1:length(all.met$GENDER)){
    gen = NA
    if(!is.na(all.met$GENDER[i])){
        gen = "male"
        if(all.met$GENDER[i] == 1){
            gen = "female"
        }
    }
    newGender = c(newGender, gen)
}
newGender = factor(newGender, levels=c("female","male")) 
all.met$GENDER = newGender



gut_set = c( "bile_duct_swab", "intest_duodenal_stent", "intest_duodenum", "intest_jejunum_swab", "panc_duct", "panc_nrml", "panc_swab", "panc_tumor")
all.met$SITEV2 = NA
all.met[all.met$SITE %in% gut_set,]$SITEV2 = "gut"
all.met[!(all.met$SITE %in% gut_set),]$SITEV2 = all.met[!(all.met$SITE %in% gut_set),]$SITE
all.met$ICDCV2 = NA
all.met[all.met$ICDC %in% c("K86.2","other"),]$ICDCV2 = "other"
all.met[!(all.met$ICDC %in% c("K86.2","other")),]$ICDCV2 = all.met[!(all.met$ICDC %in% c("K86.2","other")),]$ICDC


mouth_sites = c( "mouth_buccal","mouth_gum","mouth_saliva","mouth_tongue" )
mouth_short = c( "buccal","gum","saliva","tongue" )


sub.met = all.met[all.met$SITEV2 %in% c("gut",mouth_sites),]
sub.met$GUT = 1.0*(sub.met$SITEV2 == "gut")
sub.met$PRE_GRP = sub.met$ICDCV2

statistics.mu.pe = c()
statistics.mu.sp = c()
statistics.om.pe = c()
statistics.om.sp = c()
statistics.p0.pe = c()
statistics.p0.sp = c()
probs.mu.50 = c()
probs.mu.00 = c()
probs.om.50 = c()
probs.om.00 = c()
probs.p0.50 = c()
probs.p0.00 = c()
probs.mu.50.sp = c()
probs.mu.00.sp = c()
probs.om.50.sp = c()
probs.om.00.sp = c()
probs.p0.50.sp = c()
probs.p0.00.sp = c()
effects.mu.pe = c()
effects.mu.sp = c()
effects.om.pe = c()
effects.om.sp = c()
effects.p0.pe = c()
effects.p0.sp = c()
species = c()
muFrame = data.frame(
    MU.G1.1.E=NA, MU.G1.1.L=NA, MU.G1.1.U=NA, 
    MU.G1.2.E=NA, MU.G1.2.L=NA, MU.G1.2.U=NA, 
    MU.G1.3.E=NA, MU.G1.3.L=NA, MU.G1.3.U=NA,
    
    MU.G2.1.E=NA, MU.G2.1.L=NA, MU.G2.1.U=NA, 
    MU.G2.2.E=NA, MU.G2.2.L=NA, MU.G2.2.U=NA, 
    MU.G2.3.E=NA, MU.G2.3.L=NA, MU.G2.3.U=NA
)
omegaFrame = data.frame(
    OMEGA.G1.1.E=NA, OMEGA.G1.1.L=NA, OMEGA.G1.1.U=NA, 
    OMEGA.G1.2.E=NA, OMEGA.G1.2.L=NA, OMEGA.G1.2.U=NA, 
    OMEGA.G1.3.E=NA, OMEGA.G1.3.L=NA, OMEGA.G1.3.U=NA,
    
    OMEGA.G2.1.E=NA, OMEGA.G2.1.L=NA, OMEGA.G2.1.U=NA, 
    OMEGA.G2.2.E=NA, OMEGA.G2.2.L=NA, OMEGA.G2.2.U=NA, 
    OMEGA.G2.3.E=NA, OMEGA.G2.3.L=NA, OMEGA.G2.3.U=NA
)
p0Frame = data.frame(
    P0.G1.1.E=NA, P0.G1.1.L=NA, P0.G1.1.U=NA, 
    P0.G1.2.E=NA, P0.G1.2.L=NA, P0.G1.2.U=NA, 
    P0.G1.3.E=NA, P0.G1.3.L=NA, P0.G1.3.U=NA,
    
    P0.G2.1.E=NA, P0.G2.1.L=NA, P0.G2.1.U=NA, 
    P0.G2.2.E=NA, P0.G2.2.L=NA, P0.G2.2.U=NA, 
    P0.G2.3.E=NA, P0.G2.3.L=NA, P0.G2.3.U=NA
)

muEmpty = muFrame
omegaEmpty = omegaFrame
p0Empty = p0Frame

for(i in 1:length(sub.cnt[,1])){
    specimen = sub.cnt[i,1]
    cat(i,specimen," <---------------------------------------------------------------<\n")
    
    resp = unlist(sub.cnt[i,-1])
    sampls = colnames(sub.cnt)[-1]
    outcome = data.frame(Y = resp, SampleID = sampls)
    
    testSet = merge(outcome, sub.met, by="SampleID", all=TRUE)
    testSet = testSet[!is.na(testSet$Y),]
    
    testSet = testSet[!is.na(testSet$ICDCV2),]
    testSet = testSet[!is.na(testSet$GUT),]
    testSet = testSet[!is.na(testSet$subject_id),]
	testSet$Y = ifelse(testSet$Y<0.01,0,testSet$Y) # <1% restriction !
	
	testThisSpecies = ( sum(testSet$Y == 0) < length(testSet$Y)*0.95 )
    if(testThisSpecies && !(
			specimen %in% c( "bla" )
		)
	){ 
	
    testSet2 = testSet[,c("SampleID","Y","ICDCV2","GUT","subject_id")]
    
    testSet2$GRP = NA
    testSet2[testSet2$ICDCV2 == "C24.*" & testSet2$GUT == 0,]$GRP = 1
    testSet2[testSet2$ICDCV2 == "C25.*" & testSet2$GUT == 0,]$GRP = 2
    testSet2[testSet2$ICDCV2 == "other" & testSet2$GUT == 0,]$GRP = 3
    testSet2[testSet2$ICDCV2 == "C24.*" & testSet2$GUT == 1,]$GRP = 4
    testSet2[testSet2$ICDCV2 == "C25.*" & testSet2$GUT == 1,]$GRP = 5
    testSet2[testSet2$ICDCV2 == "other" & testSet2$GUT == 1,]$GRP = 6
	
	testSet2 = na.omit(testSet2)
    
    atab = testSet2[testSet2$Y!=0,]
    a_subject_grp_tab = data.frame( subject_id=unique(atab$subject_id), a_subject_grp=1:length(unique(atab$subject_id)) )
    atab = merge( atab, a_subject_grp_tab, by="subject_id", all=TRUE )
    
    btab = testSet2
    btab$Y = (btab$Y==0)*1
	b_subject_grp_tab = data.frame( subject_id=unique(btab$subject_id), b_subject_grp=1:length(unique(btab$subject_id)) )
	btab = merge( btab, b_subject_grp_tab, by="subject_id", all=TRUE )

    sparse_count1 = sum( as.matrix(table(atab$ICDCV2,atab$GUT)) == 0 )
    sparse_count2 = sum( as.matrix(table(btab$ICDCV2,btab$GUT)) == 0 )
    
    pr.mu.50 = pr.om.50 = pr.p0.50 = pr.mu.00 = pr.om.00 = pr.p0.00 = 
    pr.mu.50.sp = pr.om.50.sp = pr.p0.50.sp = pr.mu.00.sp = pr.om.00.sp = pr.p0.00.sp = 
    st.mu.pe = st.om.pe = st.p0.pe = st.mu.sp = st.om.sp = st.p0.sp = 
    ef.mu.pe = ef.om.pe = ef.p0.pe = ef.mu.sp = ef.om.sp = ef.p0.sp = 
    NA
    
    muEst = muEmpty
    omegaEst = omegaEmpty
    p0Est = p0Empty
    
    if(sparse_count2 == 0){
########################################################## START MODEL ABSENCE

f_mod <- function() {
    for(i in 1:N){
        Y[i] ~ dbin(p[i], 1)
        logit(p[i]) <- B[GRP[i]] + U[SBJ[i]]
    }
	
	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
	tau ~ dgamma(0.01,0.01)
    
    for(j in 1:6){
        B[j] ~ dnorm(0,0.01)
    }
    
    logit(E_P0_Oral[1]) <- B[1]
    logit(E_P0_Oral[2]) <- B[2]
    logit(E_P0_Oral[3]) <- B[3]

    logit(E_P0_Panc[1]) <- B[4]
    logit(E_P0_Panc[2]) <- B[5]
    logit(E_P0_Panc[3]) <- B[6]
    
}

dats = list(
    Y = btab$Y,
    GRP = btab$GRP, 
	SBJ = btab$b_subject_grp, 
	M = length(unique(btab$b_subject_grp)),
    N = length(btab$Y)
)

inivals <- list( list(  B=rep(0,6), tau=2 ) )

fitz <- bugs(
    dats, inits=inivals, parameters.to.save = c("E_P0_Oral","E_P0_Panc"), 
    model.file = f_mod, n.chains = 1, n.iter = MCMC_NUM+1000, n.burnin=1000, 
    digits=8 #,debug=TRUE
)

P_Oral = fitz$sims.list[["E_P0_Oral"]]
P_Panc = fitz$sims.list[["E_P0_Panc"]]

p0Mat = cbind(P_Oral, P_Panc)
p0Est = c(
    mean( na.omit(p0Mat[,1]) ), quantile( na.omit(p0Mat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,2]) ), quantile( na.omit(p0Mat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,3]) ), quantile( na.omit(p0Mat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,4]) ), quantile( na.omit(p0Mat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,5]) ), quantile( na.omit(p0Mat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,6]) ), quantile( na.omit(p0Mat[,6]), probs=c(0.025,0.975) )
)

T3 = c()
T3SP = c()
for(j in 1:dim(P_Oral)[1]){
    T3[j] = cor( x=P_Oral[j,1:3], y=P_Panc[j,1:3], method="pearson" )
    T3SP[j] = cor( x=P_Oral[j,1:3], y=P_Panc[j,1:3], method="spearman" )
}

    pr.p0.50 = mean(T3[!is.na(T3)] > 0.5)
    pr.p0.00 = mean(T3[!is.na(T3)] > 0.0)
    pr.p0.50.sp = mean(T3SP[!is.na(T3SP)] > 0.5)
    pr.p0.00.sp = mean(T3SP[!is.na(T3SP)] > 0.0)
    st.p0.pe = quantile(T3[!is.na(T3)],probs=ALPHA)
    st.p0.sp = quantile(T3SP[!is.na(T3SP)],probs=ALPHA)
    ef.p0.pe = mean(T3[!is.na(T3)])
    ef.p0.sp = mean(T3SP[!is.na(T3SP)])
    
########################################################## END MODEL ABSENCE
    }
    
    if(sparse_count1 == 0){
########################################################## START MODEL CONTINUOUS PRESENCE

f_mod <- function() {
    for(i in 1:N){
        Y[i] ~ dbeta(a[i], b[i])
        b[i] <- (1-mu[i]) * phi
        a[i] <- mu[i] * phi
        logit(mu[i]) <- B[GRP[i]] + U[SBJ[i]]
    }

    for(j in 1:6){
        B[j] ~ dnorm(0,0.01)
    }
	
	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
	tau ~ dgamma(0.01,0.01)
	
	phi <- SS*SS
	SS ~ dunif(1,100)

    logit(E_MU_Oral[1]) <- B[1]
    logit(E_MU_Oral[2]) <- B[2]
    logit(E_MU_Oral[3]) <- B[3]

    logit(E_MU_Panc[1]) <- B[4]
    logit(E_MU_Panc[2]) <- B[5]
    logit(E_MU_Panc[3]) <- B[6]
}

dats = list(
    Y = atab$Y,
    GRP = atab$GRP, 
	SIGG = atab$GUT+1,
	SBJ = atab$a_subject_grp, 
	M = length(unique(atab$a_subject_grp)),
    N = length(atab$Y)
)

inivals <- list( list(  B=rep(0,6), SS=2, tau=2 ) )

fit <- bugs(
    dats, inits=inivals, parameters.to.save = c("B","SS","E_MU_Oral","E_MU_Panc","phi"), 
    model.file = f_mod, n.chains = 1, n.iter = MCMC_NUM+1000, n.burnin=1000, 
    digits=8 #, debug=TRUE
)

M_Oral = fit$sims.list[["E_MU_Oral"]]
M_Panc = fit$sims.list[["E_MU_Panc"]]

muMat = matrix(NA,ncol=6,nrow=nrow(M_Oral))
T1 = c()
T1SP = c()
for(j in 1:dim(M_Oral)[1]){
    muMat[j,] = c( M_Oral[j,1:3]*(1-P_Oral[j,1:3]) , M_Panc[j,1:3]*(1-P_Panc[j,1:3]) )
    T1[j] = cor( x=M_Oral[j,1:3]*(1-P_Oral[j,1:3]), y=M_Panc[j,1:3]*(1-P_Panc[j,1:3]), method="pearson" )
    T1SP[j] = cor( x=M_Oral[j,1:3]*(1-P_Oral[j,1:3]), y=M_Panc[j,1:3]*(1-P_Panc[j,1:3]), method="spearman" )
}

T2 = c()
T2SP = c()
for(j in 1:dim(M_Oral)[1]){
    T2[j] = cor( x=M_Oral[j,1:3], y=M_Panc[j,1:3], method="pearson" )
    T2SP[j] = cor( x=M_Oral[j,1:3], y=M_Panc[j,1:3], method="spearman" )
}

muEst = c(
    mean( na.omit(muMat[,1]) ), quantile( na.omit(muMat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,2]) ), quantile( na.omit(muMat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,3]) ), quantile( na.omit(muMat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,4]) ), quantile( na.omit(muMat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,5]) ), quantile( na.omit(muMat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,6]) ), quantile( na.omit(muMat[,6]), probs=c(0.025,0.975) )
)

omegaMat = cbind(M_Oral, M_Panc)
omegaEst = c(
    mean( na.omit(omegaMat[,1]) ), quantile( na.omit(omegaMat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,2]) ), quantile( na.omit(omegaMat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,3]) ), quantile( na.omit(omegaMat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,4]) ), quantile( na.omit(omegaMat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,5]) ), quantile( na.omit(omegaMat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,6]) ), quantile( na.omit(omegaMat[,6]), probs=c(0.025,0.975) )
)

    pr.mu.50 = mean(T1[!is.na(T1)] > 0.5)
    pr.om.50 = mean(T2[!is.na(T2)] > 0.5)
    pr.mu.00 = mean(T1[!is.na(T1)] > 0.0)
    pr.om.00 = mean(T2[!is.na(T2)] > 0.0)
    
    pr.mu.50.sp = mean(T1SP[!is.na(T1SP)] > 0.5)
    pr.om.50.sp = mean(T2SP[!is.na(T2SP)] > 0.5)
    pr.mu.00.sp = mean(T1SP[!is.na(T1SP)] > 0.0)
    pr.om.00.sp = mean(T2SP[!is.na(T2SP)] > 0.0)
    
    st.mu.pe = quantile(T1[!is.na(T1)],probs=ALPHA)
    st.om.pe = quantile(T2[!is.na(T2)],probs=ALPHA)
    st.mu.sp = quantile(T1SP[!is.na(T1SP)],probs=ALPHA)
    st.om.sp = quantile(T2SP[!is.na(T2SP)],probs=ALPHA)
    
    ef.mu.pe = mean(T1[!is.na(T1)])
    ef.om.pe = mean(T2[!is.na(T2)])
    ef.mu.sp = mean(T1SP[!is.na(T1SP)])
    ef.om.sp = mean(T2SP[!is.na(T2SP)])
    
########################################################## END MODEL CONTINUOUS PRESENCE
    }

    if( (sparse_count1 == 0) || (sparse_count2 == 0) ){
    muFrame = rbind(muFrame, muEst)
    omegaFrame = rbind(omegaFrame, omegaEst)
    p0Frame = rbind(p0Frame, p0Est)
    
    statistics.mu.pe = c(statistics.mu.pe, st.mu.pe)
    statistics.om.pe = c(statistics.om.pe, st.om.pe)
    statistics.p0.pe = c(statistics.p0.pe, st.p0.pe)
    statistics.mu.sp = c(statistics.mu.sp, st.mu.sp)
    statistics.om.sp = c(statistics.om.sp, st.om.sp)
    statistics.p0.sp = c(statistics.p0.sp, st.p0.sp)
    effects.mu.pe = c(effects.mu.pe, ef.mu.pe)
    effects.om.pe = c(effects.om.pe, ef.om.pe)
    effects.p0.pe = c(effects.p0.pe, ef.p0.pe)
    effects.mu.sp = c(effects.mu.sp, ef.mu.sp)
    effects.om.sp = c(effects.om.sp, ef.om.sp)
    effects.p0.sp = c(effects.p0.sp, ef.p0.sp)
    probs.mu.50 = c(probs.mu.50, pr.mu.50)
    probs.om.50 = c(probs.om.50, pr.om.50)
    probs.p0.50 = c(probs.p0.50, pr.p0.50)
    probs.mu.00 = c(probs.mu.00, pr.mu.00)
    probs.om.00 = c(probs.om.00, pr.om.00)
    probs.p0.00 = c(probs.p0.00, pr.p0.00)
    probs.mu.50.sp = c(probs.mu.50.sp, pr.mu.50.sp)
    probs.om.50.sp = c(probs.om.50.sp, pr.om.50.sp)
    probs.p0.50.sp = c(probs.p0.50.sp, pr.p0.50.sp)
    probs.mu.00.sp = c(probs.mu.00.sp, pr.mu.00.sp)
    probs.om.00.sp = c(probs.om.00.sp, pr.om.00.sp)
    probs.p0.00.sp = c(probs.p0.00.sp, pr.p0.00.sp)
    species = c(species, specimen)
    }
}
}

muFrame = muFrame[-1,]
omegaFrame = omegaFrame[-1,]
p0Frame = p0Frame[-1,]

results = data.frame(
    OTU=species, 
    st.mu.pe = statistics.mu.pe, st.mu.sp = statistics.mu.sp, 
    st.om.pe = statistics.om.pe, st.om.sp = statistics.om.sp, 
    st.p0.pe = statistics.p0.pe, st.p0.sp = statistics.p0.sp, 
    ef.mu.pe = effects.mu.pe, ef.mu.sp = effects.mu.sp,
    ef.om.pe = effects.om.pe, ef.om.sp = effects.om.sp,
    ef.p0.pe = effects.p0.pe, ef.p0.sp = effects.p0.sp,
    pr.mu.pe.00 = probs.mu.00, pr.mu.sp.00 = probs.mu.00.sp, 
    pr.om.pe.00 = probs.om.00, pr.om.sp.00 = probs.om.00.sp, 
    pr.p0.pe.00 = probs.p0.00, pr.p0.sp.00 = probs.p0.00.sp, 
    pr.mu.pe.50 = probs.mu.50, pr.mu.sp.50 = probs.mu.50.sp, 
    pr.om.pe.50 = probs.om.50, pr.om.sp.50 = probs.om.50.sp, 
    pr.p0.pe.50 = probs.p0.50, pr.p0.sp.50 = probs.p0.50.sp, 
    muFrame, omegaFrame, p0Frame
)
outpath = paste0("...\\FEATURE_gut_vs_mouth_3groups_random_effects.csv")
write.csv(results, file=outpath, row.names=FALSE)











########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
#######
#######    gut vs mouth general + random effects (4 ICDC)
#######
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

options(stringsAsFactors=FALSE)
source("...\\BEINF0_utility_functions.R")
MCMC_NUM = 5000
ALPHA = 0.1

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

colsums_o = c()
for(i in 2:length(oral.cnt[1,])){
	colsums_o = c(colsums_o,sum(oral.cnt[,i]))
	oral.cnt[,i] = oral.cnt[,i]/sum(oral.cnt[,i])
}
colsumDat_o = data.frame(SampleID=colnames(oral.cnt)[-1], LIBSIZE=log(colsums_o))

colsums_t = c()
for(i in 2:length(tissue.cnt[1,])){
	colsums_t = c(colsums_t,sum(tissue.cnt[,i]))
	tissue.cnt[,i] = tissue.cnt[,i]/sum(tissue.cnt[,i])
}
colsumDat_t = data.frame(SampleID=colnames(tissue.cnt)[-1], LIBSIZE=log(colsums_t))

colsumDat = rbind(colsumDat_o,colsumDat_t)
all.met = merge(all.met,colsumDat,all=TRUE,by=c("SampleID"))



exclusionSites = c()
exclusionSites = c( exclusionSites, c("stool","stomach_swab","intest_duodenal_stent", "panc_swab") )

exclusionSet = all.met[all.met$SITE %in% exclusionSites,]$SampleID
all.met = all.met[!(all.met$SITE %in% exclusionSites),]
all.met = all.met[all.met$STUDY != "NDRI",]
all.met = all.met[!(all.met$SampleID %in% c("491762.1","491762.3","491762.4","491762.5")),]



consistent_features = intersect( oral.cnt$OTU.ID, tissue.cnt$OTU.ID )
oral.cnt = oral.cnt[oral.cnt$OTU.ID %in% consistent_features,]
tissue.cnt = tissue.cnt[tissue.cnt$OTU.ID %in% consistent_features,]



consistentSamples = unique(
	union(
		intersect(colnames(oral.cnt), all.met$SampleID),
		intersect(colnames(tissue.cnt), all.met$SampleID)
	)
)
all.met = all.met[all.met$SampleID %in% consistentSamples,]
subo.cnt = oral.cnt
subo.cnt = subo.cnt[,c("OTU.ID",colnames(subo.cnt)[colnames(subo.cnt) %in% consistentSamples])]
subt.cnt = tissue.cnt
subt.cnt = subt.cnt[,c("OTU.ID",colnames(subt.cnt)[colnames(subt.cnt) %in% consistentSamples])]

sub.cnt = merge(subt.cnt,subo.cnt,all=TRUE,by="OTU.ID")



MEDAGE = median(all.met$AGE)
MEDBMI = median(na.omit(all.met$BMI))

newGender = c()
for(i in 1:length(all.met$GENDER)){
    gen = NA
    if(!is.na(all.met$GENDER[i])){
        gen = "male"
        if(all.met$GENDER[i] == 1){
            gen = "female"
        }
    }
    newGender = c(newGender, gen)
}
newGender = factor(newGender, levels=c("female","male")) 
all.met$GENDER = newGender



gut_set = c( "bile_duct_swab", "intest_duodenal_stent", "intest_duodenum", "intest_jejunum_swab", "panc_duct", "panc_nrml", "panc_swab", "panc_tumor")
all.met$SITEV2 = NA
all.met[all.met$SITE %in% gut_set,]$SITEV2 = "gut"
all.met[!(all.met$SITE %in% gut_set),]$SITEV2 = all.met[!(all.met$SITE %in% gut_set),]$SITE
all.met$ICDCV2 = NA
all.met[all.met$ICDC %in% c("K86.2","other"),]$ICDCV2 = "other"
all.met[!(all.met$ICDC %in% c("K86.2","other")),]$ICDCV2 = all.met[!(all.met$ICDC %in% c("K86.2","other")),]$ICDC


mouth_sites = c( "mouth_buccal","mouth_gum","mouth_saliva","mouth_tongue" )
mouth_short = c( "buccal","gum","saliva","tongue" )

sub.met = all.met[all.met$SITEV2 %in% c("gut",mouth_sites),]
sub.met$GUT = 1.0*(sub.met$SITEV2 == "gut")
sub.met$PRE_GRP = sub.met$ICDC

statistics.mu.pe = c()
statistics.mu.sp = c()
statistics.om.pe = c()
statistics.om.sp = c()
statistics.p0.pe = c()
statistics.p0.sp = c()
probs.mu.50 = c()
probs.mu.00 = c()
probs.om.50 = c()
probs.om.00 = c()
probs.p0.50 = c()
probs.p0.00 = c()
probs.mu.50.sp = c()
probs.mu.00.sp = c()
probs.om.50.sp = c()
probs.om.00.sp = c()
probs.p0.50.sp = c()
probs.p0.00.sp = c()
effects.mu.pe = c()
effects.mu.sp = c()
effects.om.pe = c()
effects.om.sp = c()
effects.p0.pe = c()
effects.p0.sp = c()
species = c()
muFrame = data.frame(
    MU.G1.1.E=NA, MU.G1.1.L=NA, MU.G1.1.U=NA, 
    MU.G1.2.E=NA, MU.G1.2.L=NA, MU.G1.2.U=NA, 
    MU.G1.3.E=NA, MU.G1.3.L=NA, MU.G1.3.U=NA,
    MU.G1.4.E=NA, MU.G1.4.L=NA, MU.G1.4.U=NA,
    
    MU.G2.1.E=NA, MU.G2.1.L=NA, MU.G2.1.U=NA, 
    MU.G2.2.E=NA, MU.G2.2.L=NA, MU.G2.2.U=NA, 
    MU.G2.3.E=NA, MU.G2.3.L=NA, MU.G2.3.U=NA, 
    MU.G2.4.E=NA, MU.G2.4.L=NA, MU.G2.4.U=NA
)
omegaFrame = data.frame(
    OMEGA.G1.1.E=NA, OMEGA.G1.1.L=NA, OMEGA.G1.1.U=NA, 
    OMEGA.G1.2.E=NA, OMEGA.G1.2.L=NA, OMEGA.G1.2.U=NA, 
    OMEGA.G1.3.E=NA, OMEGA.G1.3.L=NA, OMEGA.G1.3.U=NA,
    OMEGA.G1.4.E=NA, OMEGA.G1.4.L=NA, OMEGA.G1.4.U=NA,
    
    OMEGA.G2.1.E=NA, OMEGA.G2.1.L=NA, OMEGA.G2.1.U=NA, 
    OMEGA.G2.2.E=NA, OMEGA.G2.2.L=NA, OMEGA.G2.2.U=NA, 
    OMEGA.G2.3.E=NA, OMEGA.G2.3.L=NA, OMEGA.G2.3.U=NA, 
    OMEGA.G2.4.E=NA, OMEGA.G2.4.L=NA, OMEGA.G2.4.U=NA
)
p0Frame = data.frame(
    P0.G1.1.E=NA, P0.G1.1.L=NA, P0.G1.1.U=NA, 
    P0.G1.2.E=NA, P0.G1.2.L=NA, P0.G1.2.U=NA, 
    P0.G1.3.E=NA, P0.G1.3.L=NA, P0.G1.3.U=NA,
    P0.G1.4.E=NA, P0.G1.4.L=NA, P0.G1.4.U=NA,
    
    P0.G2.1.E=NA, P0.G2.1.L=NA, P0.G2.1.U=NA, 
    P0.G2.2.E=NA, P0.G2.2.L=NA, P0.G2.2.U=NA, 
    P0.G2.3.E=NA, P0.G2.3.L=NA, P0.G2.3.U=NA, 
    P0.G2.4.E=NA, P0.G2.4.L=NA, P0.G2.4.U=NA
)

muEmpty = muFrame
omegaEmpty = omegaFrame
p0Empty = p0Frame

for(i in 1:length(sub.cnt[,1])){
    specimen = sub.cnt[i,1]
    cat(i,specimen," <---------------------------------------------------------------<\n")
    
    resp = unlist(sub.cnt[i,-1])
    sampls = colnames(sub.cnt)[-1]
    outcome = data.frame(Y = resp, SampleID = sampls)
    
    testSet = merge(outcome, sub.met, by="SampleID", all=TRUE)
    testSet = testSet[!is.na(testSet$Y),]
    
    testSet = testSet[!is.na(testSet$PRE_GRP),]
    testSet = testSet[!is.na(testSet$GUT),]
    testSet = testSet[!is.na(testSet$subject_id),]
	testSet$Y = ifelse(testSet$Y<0.01,0,testSet$Y) # IMPORTANT! <1% restriction !!!
	
	testThisSpecies = ( sum(testSet$Y == 0) < length(testSet$Y)*0.95 )
    if(testThisSpecies && !(
			specimen %in% c( "bla" )
		)
	){ 
	
    testSet2 = testSet[,c("SampleID","Y","PRE_GRP","GUT","subject_id")]
    
    testSet2$GRP = NA
    testSet2[testSet2$PRE_GRP == "C24.*" & testSet2$GUT == 0,]$GRP = 1
    testSet2[testSet2$PRE_GRP == "C25.*" & testSet2$GUT == 0,]$GRP = 2
    testSet2[testSet2$PRE_GRP == "K86.2" & testSet2$GUT == 0,]$GRP = 3
    testSet2[testSet2$PRE_GRP == "other" & testSet2$GUT == 0,]$GRP = 4
    
    testSet2[testSet2$PRE_GRP == "C24.*" & testSet2$GUT == 1,]$GRP = 5
    testSet2[testSet2$PRE_GRP == "C25.*" & testSet2$GUT == 1,]$GRP = 6
    testSet2[testSet2$PRE_GRP == "K86.2" & testSet2$GUT == 1,]$GRP = 7
    testSet2[testSet2$PRE_GRP == "other" & testSet2$GUT == 1,]$GRP = 8
    
    testSet2 = na.omit(testSet2)
    
    atab = testSet2[testSet2$Y!=0,]
    a_subject_grp_tab = data.frame( subject_id=unique(atab$subject_id), a_subject_grp=1:length(unique(atab$subject_id)) )
    atab = merge( atab, a_subject_grp_tab, by="subject_id", all=TRUE )
    
    btab = testSet2
    btab$Y = (btab$Y==0)*1
	b_subject_grp_tab = data.frame( subject_id=unique(btab$subject_id), b_subject_grp=1:length(unique(btab$subject_id)) )
	btab = merge( btab, b_subject_grp_tab, by="subject_id", all=TRUE )

    sparse_count1 = sum( as.matrix(table(atab$PRE_GRP,atab$GUT)) == 0 )
    sparse_count2 = sum( as.matrix(table(btab$PRE_GRP,btab$GUT)) == 0 )
    
    pr.mu.50 = pr.om.50 = pr.p0.50 = pr.mu.00 = pr.om.00 = pr.p0.00 = 
    pr.mu.50.sp = pr.om.50.sp = pr.p0.50.sp = pr.mu.00.sp = pr.om.00.sp = pr.p0.00.sp = 
    st.mu.pe = st.om.pe = st.p0.pe = st.mu.sp = st.om.sp = st.p0.sp = 
    ef.mu.pe = ef.om.pe = ef.p0.pe = ef.mu.sp = ef.om.sp = ef.p0.sp = 
    NA
    
    muEst = muEmpty
    omegaEst = omegaEmpty
    p0Est = p0Empty
    
    if(sparse_count2 == 0){
########################################################## START MODEL ABSENCE

f_mod <- function() {
    for(i in 1:N){
        Y[i] ~ dbin(p[i], 1)
        logit(p[i]) <- B[GRP[i]] + U[SBJ[i]]
    }
    
    for(j in 1:8){
        B[j] ~ dnorm(0,0.01)
    }
	
	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
	tau ~ dgamma(0.01,0.01)

    phi <- SS*SS
    SS ~ dunif(1,100)
    
    logit(E_P0_Oral[1]) <- B[1]
    logit(E_P0_Oral[2]) <- B[2]
    logit(E_P0_Oral[3]) <- B[3]
    logit(E_P0_Oral[4]) <- B[4]

    logit(E_P0_Panc[1]) <- B[5]
    logit(E_P0_Panc[2]) <- B[6]
    logit(E_P0_Panc[3]) <- B[7]
    logit(E_P0_Panc[4]) <- B[8]
    
}

dats = list(
    Y = btab$Y,
    GRP = btab$GRP, 
	SBJ = btab$b_subject_grp, 
	M = length(unique(btab$b_subject_grp)),
    N = length(btab$Y)
)

inivals <- list( list(  B=rep(0,8), tau=2 ) )

fitz <- bugs(
    dats, inits=inivals, parameters.to.save = c("E_P0_Oral","E_P0_Panc"), 
    model.file = f_mod, n.chains = 1, n.iter = MCMC_NUM+1000, n.burnin=1000, 
    digits=8 #,debug=TRUE
)

P_Oral = fitz$sims.list[["E_P0_Oral"]]
P_Panc = fitz$sims.list[["E_P0_Panc"]]

p0Mat = cbind(P_Oral, P_Panc)
p0Est = c(
    mean( na.omit(p0Mat[,1]) ), quantile( na.omit(p0Mat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,2]) ), quantile( na.omit(p0Mat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,3]) ), quantile( na.omit(p0Mat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,4]) ), quantile( na.omit(p0Mat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,5]) ), quantile( na.omit(p0Mat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,6]) ), quantile( na.omit(p0Mat[,6]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,7]) ), quantile( na.omit(p0Mat[,7]), probs=c(0.025,0.975) ), 
    mean( na.omit(p0Mat[,8]) ), quantile( na.omit(p0Mat[,8]), probs=c(0.025,0.975) )
)

T3 = c()
T3SP = c()
for(j in 1:dim(P_Oral)[1]){
    T3[j] = cor( x=P_Oral[j,1:4], y=P_Panc[j,1:4], method="pearson" )
    T3SP[j] = cor( x=P_Oral[j,1:4], y=P_Panc[j,1:4], method="spearman" )
}

    pr.p0.50 = mean(T3[!is.na(T3)] > 0.5)
    pr.p0.00 = mean(T3[!is.na(T3)] > 0.0)
    pr.p0.50.sp = mean(T3SP[!is.na(T3SP)] > 0.5)
    pr.p0.00.sp = mean(T3SP[!is.na(T3SP)] > 0.0)
    st.p0.pe = quantile(T3[!is.na(T3)],probs=ALPHA)
    st.p0.sp = quantile(T3SP[!is.na(T3SP)],probs=ALPHA)
    ef.p0.pe = mean(T3[!is.na(T3)])
    ef.p0.sp = mean(T3SP[!is.na(T3SP)])
    
########################################################## END MODEL ABSENCE
    }
    
    if(sparse_count1 == 0){
########################################################## START MODEL CONTINUOUS PRESENCE

f_mod <- function() {
    for(i in 1:N){
        Y[i] ~ dbeta(a[i], b[i])
        b[i] <- (1-mu[i]) * phi
        a[i] <- mu[i] * phi
        logit(mu[i]) <- B[GRP[i]] + U[SBJ[i]]
    }

	for(j in 1:M){
		U[j] ~ dnorm(0,tau)
	}
	tau ~ dgamma(0.01,0.01)
	
    for(j in 1:8){
        B[j] ~ dnorm(0,0.01)
    }

	phi <- SS*SS
	SS ~ dunif(1,100)

    logit(E_MU_Oral[1]) <- B[1]
    logit(E_MU_Oral[2]) <- B[2]
    logit(E_MU_Oral[3]) <- B[3]
    logit(E_MU_Oral[4]) <- B[4]

    logit(E_MU_Panc[1]) <- B[5]
    logit(E_MU_Panc[2]) <- B[6]
    logit(E_MU_Panc[3]) <- B[7]
    logit(E_MU_Panc[4]) <- B[8]
}

dats = list(
    Y = atab$Y,
    GRP = atab$GRP, 
	SIGG = atab$GUT+1,
	SBJ = atab$a_subject_grp, 
	M = length(unique(atab$a_subject_grp)),
    N = length(atab$Y)
)

inivals <- list( list(  B=rep(0,8), SS=2, tau=2 ) )

fit <- bugs(
    dats, inits=inivals, parameters.to.save = c("B","SS","E_MU_Oral","E_MU_Panc","phi"), 
    model.file = f_mod, n.chains = 1, n.iter = MCMC_NUM+1000, n.burnin=1000, 
    digits=8 #, debug=TRUE
)

M_Oral = fit$sims.list[["E_MU_Oral"]]
M_Panc = fit$sims.list[["E_MU_Panc"]]

muMat = matrix(NA,ncol=8,nrow=nrow(M_Oral))
T1 = c()
T1SP = c()
for(j in 1:dim(M_Oral)[1]){
    muMat[j,] = c( M_Oral[j,1:4]*(1-P_Oral[j,1:4]) , M_Panc[j,1:4]*(1-P_Panc[j,1:4]) )
    T1[j] = cor( x=M_Oral[j,1:4]*(1-P_Oral[j,1:4]), y=M_Panc[j,1:4]*(1-P_Panc[j,1:4]), method="pearson" )
    T1SP[j] = cor( x=M_Oral[j,1:4]*(1-P_Oral[j,1:4]), y=M_Panc[j,1:4]*(1-P_Panc[j,1:4]), method="spearman" )
}

T2 = c()
T2SP = c()
for(j in 1:dim(M_Oral)[1]){
    T2[j] = cor( x=M_Oral[j,1:4], y=M_Panc[j,1:4], method="pearson" )
    T2SP[j] = cor( x=M_Oral[j,1:4], y=M_Panc[j,1:4], method="spearman" )
}

muEst = c(
    mean( na.omit(muMat[,1]) ), quantile( na.omit(muMat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,2]) ), quantile( na.omit(muMat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,3]) ), quantile( na.omit(muMat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,4]) ), quantile( na.omit(muMat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,5]) ), quantile( na.omit(muMat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,6]) ), quantile( na.omit(muMat[,6]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,7]) ), quantile( na.omit(muMat[,7]), probs=c(0.025,0.975) ), 
    mean( na.omit(muMat[,8]) ), quantile( na.omit(muMat[,8]), probs=c(0.025,0.975) )
)

omegaMat = cbind(M_Oral, M_Panc)
omegaEst = c(
    mean( na.omit(omegaMat[,1]) ), quantile( na.omit(omegaMat[,1]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,2]) ), quantile( na.omit(omegaMat[,2]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,3]) ), quantile( na.omit(omegaMat[,3]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,4]) ), quantile( na.omit(omegaMat[,4]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,5]) ), quantile( na.omit(omegaMat[,5]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,6]) ), quantile( na.omit(omegaMat[,6]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,7]) ), quantile( na.omit(omegaMat[,7]), probs=c(0.025,0.975) ), 
    mean( na.omit(omegaMat[,8]) ), quantile( na.omit(omegaMat[,8]), probs=c(0.025,0.975) )
)

    pr.mu.50 = mean(T1[!is.na(T1)] > 0.5)
    pr.om.50 = mean(T2[!is.na(T2)] > 0.5)
    pr.mu.00 = mean(T1[!is.na(T1)] > 0.0)
    pr.om.00 = mean(T2[!is.na(T2)] > 0.0)
    
    pr.mu.50.sp = mean(T1SP[!is.na(T1SP)] > 0.5)
    pr.om.50.sp = mean(T2SP[!is.na(T2SP)] > 0.5)
    pr.mu.00.sp = mean(T1SP[!is.na(T1SP)] > 0.0)
    pr.om.00.sp = mean(T2SP[!is.na(T2SP)] > 0.0)
    
    st.mu.pe = quantile(T1[!is.na(T1)],probs=ALPHA)
    st.om.pe = quantile(T2[!is.na(T2)],probs=ALPHA)
    st.mu.sp = quantile(T1SP[!is.na(T1SP)],probs=ALPHA)
    st.om.sp = quantile(T2SP[!is.na(T2SP)],probs=ALPHA)
    
    ef.mu.pe = mean(T1[!is.na(T1)])
    ef.om.pe = mean(T2[!is.na(T2)])
    ef.mu.sp = mean(T1SP[!is.na(T1SP)])
    ef.om.sp = mean(T2SP[!is.na(T2SP)])
    
########################################################## END MODEL CONTINUOUS PRESENCE
    }

    if( (sparse_count1 == 0) || (sparse_count2 == 0) ){
    muFrame = rbind(muFrame, muEst)
    omegaFrame = rbind(omegaFrame, omegaEst)
    p0Frame = rbind(p0Frame, p0Est)
    
    statistics.mu.pe = c(statistics.mu.pe, st.mu.pe)
    statistics.om.pe = c(statistics.om.pe, st.om.pe)
    statistics.p0.pe = c(statistics.p0.pe, st.p0.pe)
    statistics.mu.sp = c(statistics.mu.sp, st.mu.sp)
    statistics.om.sp = c(statistics.om.sp, st.om.sp)
    statistics.p0.sp = c(statistics.p0.sp, st.p0.sp)
    effects.mu.pe = c(effects.mu.pe, ef.mu.pe)
    effects.om.pe = c(effects.om.pe, ef.om.pe)
    effects.p0.pe = c(effects.p0.pe, ef.p0.pe)
    effects.mu.sp = c(effects.mu.sp, ef.mu.sp)
    effects.om.sp = c(effects.om.sp, ef.om.sp)
    effects.p0.sp = c(effects.p0.sp, ef.p0.sp)
    probs.mu.50 = c(probs.mu.50, pr.mu.50)
    probs.om.50 = c(probs.om.50, pr.om.50)
    probs.p0.50 = c(probs.p0.50, pr.p0.50)
    probs.mu.00 = c(probs.mu.00, pr.mu.00)
    probs.om.00 = c(probs.om.00, pr.om.00)
    probs.p0.00 = c(probs.p0.00, pr.p0.00)
    probs.mu.50.sp = c(probs.mu.50.sp, pr.mu.50.sp)
    probs.om.50.sp = c(probs.om.50.sp, pr.om.50.sp)
    probs.p0.50.sp = c(probs.p0.50.sp, pr.p0.50.sp)
    probs.mu.00.sp = c(probs.mu.00.sp, pr.mu.00.sp)
    probs.om.00.sp = c(probs.om.00.sp, pr.om.00.sp)
    probs.p0.00.sp = c(probs.p0.00.sp, pr.p0.00.sp)
    species = c(species, specimen)
    }
}
}

muFrame = muFrame[-1,]
omegaFrame = omegaFrame[-1,]
p0Frame = p0Frame[-1,]

results = data.frame(
    OTU=species, 
    st.mu.pe = statistics.mu.pe, st.mu.sp = statistics.mu.sp, 
    st.om.pe = statistics.om.pe, st.om.sp = statistics.om.sp, 
    st.p0.pe = statistics.p0.pe, st.p0.sp = statistics.p0.sp, 
    ef.mu.pe = effects.mu.pe, ef.mu.sp = effects.mu.sp,
    ef.om.pe = effects.om.pe, ef.om.sp = effects.om.sp,
    ef.p0.pe = effects.p0.pe, ef.p0.sp = effects.p0.sp,
    pr.mu.pe.00 = probs.mu.00, pr.mu.sp.00 = probs.mu.00.sp, 
    pr.om.pe.00 = probs.om.00, pr.om.sp.00 = probs.om.00.sp, 
    pr.p0.pe.00 = probs.p0.00, pr.p0.sp.00 = probs.p0.00.sp, 
    pr.mu.pe.50 = probs.mu.50, pr.mu.sp.50 = probs.mu.50.sp, 
    pr.om.pe.50 = probs.om.50, pr.om.sp.50 = probs.om.50.sp, 
    pr.p0.pe.50 = probs.p0.50, pr.p0.sp.50 = probs.p0.50.sp, 
    muFrame, omegaFrame, p0Frame
)
outpath = paste0("...\\FEATURE_gut_vs_mouth_4groups_random_effects.csv")
write.csv(results, file=outpath, row.names=FALSE)
