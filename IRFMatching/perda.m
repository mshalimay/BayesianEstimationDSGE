function [wlf,IRt] = perda(param,set,IRFempirica,W,mod,choque,tun)
%param sao os parametros sendo estimados (vetor)
%set sao os parametros mantidos fixos (vetor)
%IRFempirica eh uma matriz (nvar x T) contendo as respostas das nvar variaveis
% W eh a inversa da matriz diagonal (nvar x T x nvar x T) dos desvios padroes de cada impulso do vetor IRFempirica na diagonal principal
%mod eh o arquivo com o modelo original (so para facilitar indexacao)
%choque eh o choque exogeno que esta sendo usado em formato string com o mesmo nome que no arquivo mod

%%*************************************************************************
%Solucao do modelo (GENSYS)
%**************************************************************************
[g0 g1 PSI PI]=model_prog(param,set);
neq=size(g0,1);
C=zeros(neq,1);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);
%%
%*************************************************************************
%Calculo da funcao de resposta ao impulso
%*************************************************************************
T=size(IRFempirica,1);
IRt=zeros(T,neq);
%Ajustando os choques para coincidir com as dinamicas do paper
tunning=zeros(size(impact,1),1);
tunning(find(mod.eps=='epsvet'))=1.25;
tunning(find(mod.eps=='epsrstar'))=-5.5661;
tunning(find(mod.eps=='epscstar'))=2.4;
tunning(find(mod.eps=='epsmm'))=tun;
%Calcula a impulso resposta
[IR]=ir_gensys(G1,impact,T,tunning); %tunning opcional para ajustar o tamanho do choque na IRF. Se vazio, tunning=1
%transforma Impulsos de inflacao e juros em anuais (como no paper)
IR(:,find(mod.Y=='pidd'),:)=4*IR(:,find(mod.Y=='pidd'),:);
IR(:,find(mod.Y=='pif'),:)=4*IR(:,find(mod.Y=='pif'),:);
IR(1,find(mod.Y=='pidd'),:)=0; %robada do SW
%coleta as impulso resposta relevantes
IRt=IR(1:T,[find(mod.Y=='pidd'),find(mod.Y=='pif')],find(mod.eps==choque));
%%
%*************************************************************************
%Funcao perda
%*************************************************************************
%plot([IRt(:,1),IRFempirica(:,1)])
%plot([IRt(:,2),IRFempirica(:,2)])

desv = IRFempirica(:) - IRt(:);

wlf = 1000*desv'*W*desv;
