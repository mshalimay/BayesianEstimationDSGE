boot_i<-function(presmply,xboot,u,lagsy,inicio,BB,Bx,nIRF){

uboot=as.data.frame(u-colMeans(u))
uboot=sapply(uboot,function(x) sample(x,replace = TRUE,size = (nrow(u)+lagsy)))

yboot=matrix(0,nrow(uboot),ncol(u))

yboot[1:lagsy,]=presmply

for(i in (lagsy+1):nrow(uboot)){
  yboot[i,]=BB%*%c(t(yboot[(i-1):(i-lagsy),]))+Bx%*%xboot[i,]+uboot[i,]
}

colnames(yboot)=colnames(presmply)
#*************************************************
#*Var na FR
#*************************************************
varred=rfvar3(yboot,lagsy,xboot,FALSE,NULL,NULL,NULL,NULL)
B_boot=varred$By
uu=varred$u
sigmaa=var(uu)

Bx_boot=varred$Bx
seq=seq(from = 0,by = ncol(Bx_boot),to=ncol(Bx_boot)*nrow(Bx_boot))
Bx_boot=Reduce('rbind',lapply(1:nrow(Bx_boot),function(i) Bx_boot[(seq[i]+1):(seq[i+1])]))
rownames(Bx_boot)=rownames(B_boot)


#*************************************************
#IDENTIFICACAO DO SVAR
##OBS.: Seja e_t erro da FR e u_t erro estru. . Autores usam RATS, que faz: G u_t = e_t
## Eviews faz: A e_t = B u_t . Codigos aqui fazem: A0^(-1) u_t = e_t 
#(Ou seja, saindo de um VAR escrito como A0 y_t = A1 yt-1 + ...+ Ap yt-p + u_t )
# Identificacao como os autores, impondo restricoes sobre a matriz G. Depois transformo G em A0 para
#computar a verossimilhanca do SVAR.
idmat=matrix(FALSE,ncol(sigmaa),ncol(sigmaa))
rownames(idmat)=rownames(sigmaa)
colnames(idmat)=colnames(sigmaa)

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
H0=numHess(SVARlh0,pvec,idmat=idmat,sigma=sigmaa,T=nrow(uu))
svar=csminwelNew(SVARlh0,x0=pvec,H0=H0,idmat=idmat,sigma=sigmaa,T=nrow(uu),crit = 1e-20,nit=10000)

#Constroi G e A0 = G^(-1)
n <- dim(idmat)[1]
G_boot <- matrix(0,n,n)
G_boot[idmat] <- svar$xh
diag(G_boot) <- 1 #efeitos unitarios da variavel em si mesma
rownames(G_boot)=rownames(sigmaa)
colnames(G_boot)=colnames(sigmaa)
G_boot['CPI','REX']=0; #robadinha SW
#A0=solve(G_boot)

#Structural decomposition
IRF_boot=impulsdtrf(B_boot,smat=G_boot,nIRF)

return(list(B_boot=B_boot,Bx_boot=Bx_boot,G_boot=G_boot,IRF_boot=IRF_boot))
}
