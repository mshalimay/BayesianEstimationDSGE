rm(list=ls())
library(ggplot2)
library(vars)
require(zoo)
library(multicore)
library(boot)
library(parallel)
library(foreach)
library(doParallel)
library(reshape)
library(R.matlab)


cl <- makeCluster(3)
registerDoParallel(cl,cores=3)
options(scipen=999)

setwd("C:/Users/Moises/Desktop/papers para replicar/Estimacao/R")
source("rfvar3.r");source("impulsdtrf.r");
source("SVARlh.r");source("SVARlh0.r")
source("csminwelNew.r");source("numgrad.r");source("csminit.r");source("numHess.r");source("bfgsi.r")
source("boot_i.r")

#*************************************************
#*Organizando matrizes do VAR
#*************************************************
dados<-read.table("input.csv",header=TRUE,sep=",",row.names=1)

x=dados[c("YUS","CPIUS","RUS","WP")]
y=setdiff(dados,x)

lagsy = 2
lagsx = 3
##onde comeca a amostra (-2 pq o programa do SIMS ja tira duas obs.)
iniciosims=which(row.names(y)=='77Q1')-lagsy
inicio=iniciosims+lagsy


##CRIA AS VARIAVEIS COM LAG
xx=Reduce(cbind,sapply(0:lagsx,function(i) lag(zoo(x),-i)))
yy=Reduce(cbind,sapply(1:lagsy,function(i) lag(zoo(y),-i)))

#CRIA OS FATORES DE SAZONALIDADE LEMBRAR DE AJUSTAR CONFORME A PLANILHA DE INPUTS PARA QUE CAIAM NO TRIMESTRE CERTO
seas1=rep(c(1,0,0,0),nrow(dados)/4); seas2=rep(c(0,1,0,0),nrow(dados)/4); seas3=rep(c(0,0,1,0),nrow(dados)/4)
seas=cbind(seas1,seas2,seas3)

#matriz de regressores
X=as.matrix(cbind(seas,xx)) #se for usar o prog. do sims
X=X[iniciosims:nrow(X),]

XX=cbind(1,as.matrix(cbind(seas,xx,yy))) #na mao ou por LM
XX=XX[inicio:nrow(XX),]
#matriz de explicadas
Y=as.matrix(y[iniciosims:nrow(y),])
YY=as.matrix(y[inicio:nrow(y),])

#BETA=solve(crossprod(X,X),crossprod(X,Y))# <<---CALCULO NA MAO
#varred=lm(YY~0+XX) #calculo por OLS

#*************************************************
#*Var na FR
#*************************************************
varred=rfvar3(Y,2,X,TRUE,NULL,NULL,NULL,NULL)
B=varred$By
u=varred$u

Bx=varred$Bx
seq=seq(from = 0,by = ncol(Bx),to=ncol(Bx)*nrow(Bx))
Bx=Reduce('rbind',lapply(1:nrow(Bx),function(i) Bx[(seq[i]+1):(seq[i+1])]))
rownames(Bx)=rownames(B)

sigma=var(u)

#*************************************************
#IDENTIFICACAO DO SVAR
##OBS.: Seja e_t erro da FR e u_t erro estru. . Autores usam RATS, que faz: G u_t = e_t
## Eviews faz: A e_t = B u_t . Codigos aqui fazem: A0^(-1) u_t = e_t 
#(Ou seja, saindo de um VAR escrito como A0 y_t = A1 yt-1 + ...+ Ap yt-p + u_t )
# Identificacao como os autores, impondo restricoes sobre a matriz G. Depois transformo G em A0 para
#computar a verossimilhanca do SVAR.

idmat=matrix(FALSE,ncol(sigma),ncol(sigma))
rownames(idmat)=rownames(sigma)
colnames(idmat)=colnames(sigma)

idmat['Y','Y']=TRUE; #PRODUTO NAO RECEBE EFEITOS CONTEMPORANEOS
idmat['XMY','XMY']=TRUE; idmat['XMY','Y']=TRUE; #EXPORTACOES SO SAO AFETADAS PELO PRODUTO
idmat['CPI','CPI']=TRUE;idmat['CPI','XMY']=TRUE; idmat['CPI','Y']=TRUE;idmat['CPI','REX']=TRUE; #cambio, produto e exportacoes afetam inflacao domestica
idmat['R','CPI']=TRUE;idmat['R','R']=TRUE #SO INFLACAO DOMESTICA AFETA JUROS
idmat['REX',]=TRUE; idmat['REX','CPIM']=FALSE; 
idmat['CPIM',]=TRUE #inflacao de importados sofre efeitos de todos

diag(idmat)<-FALSE #elementos da diagonal nao precisam de identificacao
##************************************************

#PARAMETROS E CHUTE INICIAL***********************
pvec=idmat
pvec=idmat*0 #chute do SW
diag(pvec)<-1 #G COM DIAG = 1
pvec['R','CPI']=0.05;pvec['CPI','REX']=-0.25; ##chute do SW
pvec=pvec[idmat==TRUE]
#*************************************************

#Otimizacao
  H0=numHess(SVARlh0,pvec,idmat=idmat,sigma=sigma,T=nrow(u))
svar=csminwelNew(SVARlh0,x0=pvec,H0=H0,idmat=idmat,sigma=sigma,T=nrow(u),crit = 1e-40,nit=10000)
svarr=optim(pvec,SVARlh0,idmat=idmat,sigma=sigma,T=nrow(u),control=list(maxit=100000,reltol=1e-40))

#Constroi G e A0 = G^(-1)
n <- dim(idmat)[1]
G <- matrix(0,n,n)
G[idmat] <- svar$xh
diag(G) <- 1 #efeitos unitarios da variavel em si mesma
rownames(G)=rownames(sigma)
colnames(G)=colnames(sigma)
G['CPI','REX']=0; #robadinha SW
#A0=solve(G)

#Structural decomposition
IRF=impulsdtrf(B,smat=G,20)

IRF_r_cpi=IRF[,'R',]['CPI',]
IRF_r_cpim=IRF[,'R',]['CPIM',]
IRF_rex_cpi=IRF[,'REX',]['CPI',]

IRF_rex_cpim=IRF[,'REX',]['CPIM',]

plot(1:20,IRF_rex_cpi)
abline(h=0)

#******************************************
#Bootstrap
#******************************************
##CRIA AS VARIAVEIS COM LAG
xx=Reduce(cbind,sapply(0:lagsx,function(i) lag(zoo(x),-i)))
yy=Reduce(cbind,sapply(1:lagsy,function(i) lag(zoo(y),-i)))

#CRIA OS FATORES DE SAZONALIDADE LEMBRAR DE AJUSTAR CONFORME A PLANILHA DE INPUTS PARA QUE CAIAM NO TRIMESTRE CERTO
seas1=rep(c(1,0,0,0),nrow(dados)/4); seas2=rep(c(0,1,0,0),nrow(dados)/4); seas3=rep(c(0,0,1,0),nrow(dados)/4)
seas=cbind(seas1,seas2,seas3)

#matriz de exogenas
xboot=cbind(as.matrix(cbind(seas,xx)),1)
xboot=xboot[iniciosims:nrow(xboot),]

#condicoes iniciais dadas
presmply=as.matrix(y[iniciosims:(iniciosims+lagsy-1),])

#reescreve a matriz B das endogenas do VAR no formato adequado
BB=B[,,1]
for(i in 2:dim(B)[[3]]){
  BB=cbind(BB,B[,,i])  
}


#bootstraping: 1000 replicacoes
bootstrap=list()
while(length(bootstrap)<1000){
  end=1000-length(bootstrap)
  print(length(bootstrap))
  boot<-foreach(i=1:end,.errorhandling = "remove") %dopar% {
  boot_i(presmply,xboot,u,lagsy,iniciosims,BB,Bx,20)}
  bootstrap=c(bootstrap,boot)
}
#Bootstraping: corrige o vies pela estimativa do bootstrap anterior e simula 2000 vezes

#corrigind vies
b=Reduce('+',lapply(bootstrap,function(x) x$B_boot))/length(bootstrap)
bx=Reduce('+',lapply(bootstrap,function(x) x$Bx_boot))/length(bootstrap)
viesb=b-B
viesbx=bx-Bx
Bstar=B-viesb
Bxstar=Bx-viesbx

BBstar=Bstar[,,1]
for(i in 2:dim(Bstar)[[3]]){
  BBstar=cbind(BBstar,Bstar[,,i])  
}
#corrige os residuos
xstar=xboot[(lagsy+1):nrow(xboot),]
ystar=as.matrix(yy[(inicio-1):nrow(yy),])
ustar=as.matrix(y[inicio:nrow(y),]-ystar%*%t(BBstar)-xstar%*%t(Bxstar))

bootstrap=list()
# 2000 replicacoes
while(length(bootstrap)<2000){
  end=2000-length(bootstrap)
  print(length(bootstrap))
  boot<-foreach(i=1:end,.errorhandling = "remove") %dopar% {
    boot_i(presmply,xboot,ustar,lagsy,iniciosims,BBstar,Bxstar,20)}
  bootstrap=c(bootstrap,boot)
}
#salva o exercicio de simulacao
save(bootstrap,file='bootstrap')

#******************************************
#Plot e desvio padrao das IRF's
#******************************************
#load nos resultados do bootstrap
load('bootstrap')
tunningr = 0.35 ##para fittar com o tamanho do choque monetario em SW
tunningrex=1.7  ##para fittar com o tamanho do choque cambial em SW

#IR do choque monetario
IRFr=lapply(bootstrap,function(x) t(x$IRF[,'R',])*tunningr) ##0.35 PARA FITTAR COM O TAMANHO DO CHOQUE QUE SW DaO
IRFrbar=Reduce('+',IRFr)/length(IRFr)
IRFrsd=sqrt(Reduce('+',lapply(IRFr,function(x) (x-IRFrbar)^2))/length(IRFr))

#IRF do choque de cambio
IRFrex=lapply(bootstrap,function(x) t(x$IRF[,'REX',])*tunningrex)
IRFrexbar=Reduce('+',IRFrex)/length(IRFrex)
IRFrexsd=sqrt(Reduce('+',lapply(IRFrex,function(x) (x-IRFrexbar)^2))/length(IRFrex))

#Intervalo de confianca da IRF do choque monetario
IRFr=t(IRF[,'R',])*tunningr #0.35 eh para fittar com o choque do SW
upperR=IRFr+2*IRFrsd
lowerR=IRFr-2*IRFrsd

#Intervalo de confianca do choque de cambio
IRFrex=t(IRF[,'REX',])*tunningrex
upperrex=IRFrex+2*IRFrexsd
lowerrex=IRFrex-2*IRFrexsd

#Exportando para matlab
IRFrcpi=IRFr[,c('CPI','CPIM')]
IRFrcpisd=IRFrsd[,c('CPI','CPIM')]
IRFrexcpi=IRFrex[,c('CPI','CPIM')]
IRFrexcpisd=IRFrsd[,c('CPI','CPIM')]

filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRFr", ".mat", sep="")
writeMat(filename, IRFr=IRFrcpi)
filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRFrsd", ".mat", sep="")
writeMat(filename, IRFrsd=IRFrcpisd)

filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRFrex", ".mat", sep="")
writeMat(filename, IRFrex=IRFrexcpi)
filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRFrexsd", ".mat", sep="")
writeMat(filename, IRFrexsd=IRFrexcpisd)

filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRr", ".mat", sep="")
writeMat(filename, IRr=IRFr[,c('R','REX','Y','XMY','CPI','CPIM')])
filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRrsd", ".mat", sep="")
writeMat(filename, IRrsd=IRFrsd[,c('R','REX','Y','XMY','CPI','CPIM')])

filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRrex", ".mat", sep="")
writeMat(filename, IRrex=IRFrex[,c('R','REX','Y','XMY','CPI','CPIM')])
filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRrexsd", ".mat", sep="")
writeMat(filename, IRrexsd=IRFrexsd[,c('R','REX','Y','XMY','CPI','CPIM')])


filename <- paste("C:/Users/Moises/Desktop/papers para replicar/Estimacao/Matlab/IRFrsd", ".mat", sep="")
writeMat(filename, IRFrsd=IRFrcpisd)


#Plotando************************************************************
df<-data.frame(t=1:20,irf=IRFrex[,'CPIM'],u=upperrex[,'CPIM'],d=lowerrex[,'CPIM'])
df<-melt(df,id="t")

graph=ggplot(df,aes(x=t,y=value,color=variable))+geom_line(aes(linetype=variable))+
  scale_color_manual(values = c("dodgerblue2", "red", "red")) +  
  scale_linetype_manual(values=c("solid","dotted","dotted"))

graph+theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none")
#*********************************************************************




