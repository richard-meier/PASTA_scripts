
library(gamlss)
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
  ##print("Utilizing hacked openbugs script to increase chain precision!")
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

options(stringsAsFactors=FALSE)

dbeta2=function(x,mu,phi){
	alp = mu*phi
	bet = (1-mu)*phi
	dbeta(x,alp,bet)
}

p.logit = function(x){ log((1+x)/(1-x)) }
p.sigmoid = function(y){ (exp(y)-1)/(exp(y)+1) }





startTime = Sys.time()
setwd("...\\new_sims")

TCONC = 0.0
terr = 0.001
R = 1000
MCMC_NUM = 5000
MCMC_BURN = 500

EACH_NUMS = c(20,40,60,80)

for( eee in 1:length(EACH_NUMS)){

cat("\t",eee," NS ======================== \n")
each_num = EACH_NUMS[eee]

DELTAS = c(0.8,1.6,2.4,3.2)
for( ddd in 1:length(DELTAS) ){

set.seed(123)
cat("\t\t",ddd," DELTA ======================== \n")
delta = DELTAS[ddd]

c_tconc = p.sigmoid( p.logit(TCONC) + delta )

p0A = NA
p0B = NA
mmu = NA
phi = NA
YS = NA

fitps = list()
fitms = list()
true.P0 = matrix(ncol=8,nrow=R)
pred.P0 = matrix(ncol=8,nrow=R)
QMP = matrix(ncol=99,nrow=R)
QMS = matrix(ncol=99,nrow=R)

for(k in 1:R){

cat(k,"\t")

rerun = TRUE

while(rerun){
p0A = rbeta(4, 7.35, 0.49)
p0B = rbeta(4, 7.35, 0.49)

true_t = cor(p0A,p0B,method="pearson")
while( abs(c_tconc - true_t) > terr ){
	p0A = rbeta(4, 7.35, 0.49)
	p0B = rbeta(4, 7.35, 0.49)
	true_t = cor(p0A,p0B,method="pearson")
}
P0S = rep(c(p0A,p0B),each=each_num)
N = length(P0S)

YS = rbinom(n=N,size=1,prob=P0S)
GRP = rep(1:8,each=each_num)
btab = data.frame(Y=YS,GRP=GRP,p0.true=P0S)
sparse_count = sum(table(factor(btab$Y,levels=c(0,1)),btab$GRP) == 0)

limit = 1
while(sparse_count > 3 & limit <= 1000){
	limit = limit + 1
	YS = rbinom(n=N,size=1,prob=P0S)
	GRP = rep(1:8,each=each_num)
	btab = data.frame(Y=YS,GRP=GRP,p0.true=P0S)
	sparse_count = sum(table(factor(btab$Y,levels=c(0,1)),btab$GRP) == 0)
}

if(sparse_count <= 3){ rerun=FALSE }
}

f_mod <- function() {
	for(i in 1:N){
		Y[i] ~ dbin(p[i], 1)
		logit(p[i]) <- B[GRP[i]]
	}

	for(j in 1:8){
		B[j] ~ dnorm(0,0.01)
	}
	
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
	N = length(btab$Y)
)

inivals <- list( list(  B=rep(0,8), SS=2 ) )

fitz <- bugs(
	dats, inits=inivals, parameters.to.save = c("B","E_P0_Oral","E_P0_Panc"), 
	model.file = f_mod, n.chains = 1, n.iter = MCMC_NUM+MCMC_BURN, n.burnin=MCMC_BURN, 
	digits=8 #, debug=TRUE
)

smt = fitz$summary
predP0 = smt[substr(row.names(smt),1,4) == "E_P0",1]

P_Oral = fitz$sims.list[["E_P0_Oral"]]
P_Panc = fitz$sims.list[["E_P0_Panc"]]

TP = c()
TS = c()
for(i in 1:dim(P_Oral)[1]){
	TP[i] = cor( x=P_Oral[i,1:4], y=P_Panc[i,1:4], method="pearson" )
	TS[i] = cor( x=P_Oral[i,1:4], y=P_Panc[i,1:4], method="spearman" )
}
qvecp = quantile(TP, probs = seq(from=0.01,to=0.99,by=0.01))
qvecs = quantile(TS, probs = seq(from=0.01,to=0.99,by=0.01))
QMP[k,] = qvecp
QMS[k,] = qvecs

true.P0[k,] = c(p0A,p0B)
pred.P0[k,] = predP0
print(Sys.time() - startTime)

}

colnames(QMP) = paste("q",seq(from=0.01,to=0.99,by=0.01),sep=".")
sim = list(QMP=QMP, QMS=QMS, P0.TRUE=true.P0, P0.PRED=pred.P0)

outpath = paste0("sim_bernP0_ft_","n",each_num,"_T",TCONC*100,"_delta",delta*10,"_grps4_spct3_pe.Rdata")
save(sim, c_tconc, delta, file=outpath)

}

}

