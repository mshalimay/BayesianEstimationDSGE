clear
close all
tun = 2 %para rhom=0.7
%tun=5.5 %para rhom=0.55
%**********************************************************
% Resolve o modelo analiticamente
%**********************************************************

%parametros
[param0,set] = parameters_est;

%Define o modelo e calcula as matrizes do gensys analiticamente
mod = modelsims(param0,set);
%dimensao do modelo
neq=length(mod.f);
%Gerar uma funcao que retorna as matrizes do Gensys 
%essa linha criara um file ".m" com as matrizes simbolicas que
%serao usadas na solucao do modelo (nome do arquivo criado: model_prog)
%aumentar eficiencia: nao precisa computar as derivadas em toda iteracao
model_func_sims(mod);

%Valores iniciais dos parametros
param0 = struct2array(param0);
set = struct2array(set);

%*************************************************************************
%Load das IRFs e calculo da matriz de ponderacao
%*************************************************************************

%Primeira coluna eh CPI domehstico e segunda coluna o estrangeiro
%empiricas
load('IRFr'); load('IRFrsd'); load('IRFrex'); load('IRFrexsd');
T=100;
%Matrizes de ponderacao
IRFrvar=IRFrsd.^2;
IRFrvar=IRFrvar(:);
Wr=diag(IRFrvar,0);

IRFrexvar=IRFrexsd.^2;
IRFrexvar=IRFrexvar(:);
Wrex=diag(IRFrexvar,0);


%%*************************************************************************
%%Calcula as IRF's teoricas
%%*************************************************************************
Matrizes do Gensys com valores iniciais
param0=[0.9018,0.8806,0,0.5365];
[g0 g1 PSI PI]=model_prog(param0,set);
C=zeros(neq,1);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);

%Calculo da funcao de resposta ao impulso
T=30;
IRv=zeros(T,neq);IRrstar=zeros(T,neq);IRcstar=zeros(T,neq);IRm=zeros(T,neq);
%Ajustando os choques para coincidir com as dinamicas do paper
%tunning=ones(4,1);
tunning=zeros(size(impact,1),1);
tunning(find(mod.eps=='epsvet'))=1.25;
tunning(find(mod.eps=='epsrstar'))=-5.5661;
tunning(find(mod.eps=='epscstar'))=2.4;
tunning(find(mod.eps=='epsmm'))=tun;
%Calcula a impulso resposta
[IR]=ir_gensys(G1,impact,T,tunning) ;%tunning opcional para ajustar o tamanho do choque na IRF. Se vazio, tunning=1
%transforma Impulsos de inflacao e juros em anuais (como no paper)
IR(:,find(mod.Y=='pidd'),:)=4*IR(:,find(mod.Y=='pidd'),:);
IR(:,find(mod.Y=='pif'),:)=4*IR(:,find(mod.Y=='pif'),:);
IR(:,find(mod.Y=='pi'),:)=4*IR(:,find(mod.Y=='pi'),:);
IR(:,find(mod.Y=='rr'),:)=4*IR(:,find(mod.Y=='rr'),:);
%coleta impulsos de cada choque em vetores especificos
IRrstar=IR(:,:,find(mod.eps=='epsrstar'));
IRv=IR(:,:,find(mod.eps=='epsvet'));
IRcstar=IR(:,:,find(mod.eps=='epscstar'));
IRm=IR(:,:,find(mod.eps=='epsmm'));

%IRFs usadas na otimizacao
IRFrt=IRm(1:T,[find(mod.Y=='pidd'),find(mod.Y=='pif')])
IRFrext=IRrstar(1:T,[find(mod.Y=='pidd'),find(mod.Y=='pif')]);

%% 
%**************************************************************************
%******************* IRF MATCHING******************************************
%**************************************************************************

%***************************************
%Funcao perda da IR ao choque monetario
%***************************************

%options da otimizacao
options = optimset('fmincon'); 
options = optimset(options, 'TolFun', 1e-40);

wr=pinv(Wr);

param0=[0.9, 0.9, 0, 0.45];
perda(param0,set,IRFr,wr,mod,'epsmm',tun);
objh = @(param) perda(param,set,IRFr,wr,mod,'epsmm',tun);
[param_r,wlf_r,exitflag,output,lambda,grad,hessian_r] = fmincon(objh, param0,[],[],[],[],[.01,.01,0,0],[.99,.99,.99,.99],[],options);
dp_r=sqrt(diag(hessian_r^(-1)));
prob_r=chi2cdf(wlf_r/1000,19);
[lixo,IRtr]=perda(param_r,set,IRFr,wr,mod,'epsmm',tun)

%--------------------------------------------------------------------------
%Esse pedaco testa condicoes iniciais para a minimizacao da funcao
%perda quando o choque eh o choque monetario.

estimativas=zeros(980,5);
wr=pinv(Wr)
matlabpool ('open');
parfor i=1:980;
	 param0=[i/1000, i/1000, (981-i)/1000, (981-i)/1000];
	 objh = @(param) perda(param,set,IRFr,wr,mod,'epsmm',tun);
	 [param_r,wlf_r] = fmincon(objh, param0,[],[],[],[],[.01,.01,0,0],[.99,.99,.99,.99],[],options);
	 estimativas(i,:)=[param_r wlf_r/1000];
 end
matlabpool close;
save('estimativas',estimativas);
%--------------------------------------------------------------------------

%***************************************
%Funcao perda da IR ao choque de cambio
%***************************************
param0=[0.9,0.9,0.5,0.5];
objh = @(param) perda(param,set,IRFrex,pinv(Wrex),mod,'epsrstar',0);
[param_rex,wlf_rex,exitflag,output,lambda,grad,hessian_rex]  = fmincon(objh, param0,[],[],[],[],[.01,.01,0,0],[.99,.99,.99,.99],[],options);
dp_rex=diag(hessian_rex^(-1));
prob_rex=chi2cdf(wlf_rex/1000,19);
[lixo,IRtrex]=perda(param_rex,set,IRFrex,pinv(Wrex),mod,'epsrstar',0);
%%
save('IR_teorica.mat','IRtr','IRtrex')

%--------------------------------------------------------------------------
%Esse pedaco testa condicoes iniciais para a minimizacao da funcao
%perda quando o choque eh o choque de cambio.

estimativas=zeros(980,5);
wr=pinv(Wrex)
matlabpool ('open');
parfor i=1:980;
 param0=[i/1000, i/1000, (981-i)/1000, (981-i)/1000];
 objh = @(param) perda(param,set,IRFrex,Wrex,mod,'epsrstar',tun);
 [param_rex,wlf_rex] = fmincon(objh, param0,[],[],[],[],[.01,.01,0,0],[.99,.99,.99,.99],[],options);
 estimativas(i,:)=[param_rex wlf_rex/1000];
 end
matlabpool close;
save('estimativas',estimativas);
%--------------------------------------------------------------------------


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Essa secao foi feita pra tentar fitar com as IRFs do paper (terceira linha fig1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************************************************************
%  %Esse pedaco foi feito para tentar fittar com a IRF do choque de cambio 
%  no paper graficamente
%**************************************************************************
%  inicio=1
%  x=IRFrext(inicio:inicio+20,1)
%  x(1,1)=0
%  um=IRFrex+2*IRFrexsd
%  lm=IRFrex-2*IRFrexsd
%  figure
%  plot([IRFrex(:,1),um(:,1),lm(:,1),x(1:20)])
% 
% inicio=1
% x=IRFrext(inicio:inicio+20,2)
% plot(x)
% um=IRFrex+2*IRFrexsd
% lm=IRFrex-2*IRFrexsd
% figure
% plot([IRFrex(:,2),um(:,2),lm(:,2),x(1:20)])

%%*************************************************************************
 %Esse pedaco foi feito para tentar fittar com a IRF do choque de mon. 
 %%no paper graficamente
%%************************************************************************** 
%inicio=1
% x=IRFrt(inicio:inicio+20,1)
% x(1,1)=0
% um=IRFr+2*IRFrsd
% lm=IRFr-2*IRFrsd
% figure
% plot([IRFr(:,1),um(:,1),lm(:,1),x(1:20)])
%
%inicio=1
%x=IRFrt(inicio:inicio+20,2)
%plot(x)
%um=IRFr+2*IRFrsd
%lm=IRFr-2*IRFrsd
%figure
%plot([IRFr(:,2),um(:,2),lm(:,2),x(1:20)])
%

%**************************************************************************
%Funcao de ajuste para tentar fittar com as IRFs do paper, minimizando a
%funcao perda em relacao ao tamanho do choque e coeficiente autoregressivo
%%sob os parametros 'verdadeiros' estimados por SW
%**************************************************************************
%options = optimset('fmincon'); 
%options = optimset(options, 'TolFun', 1e-20, 'display', 'iter');
%fittando as IRFs com as do SW
%para o choque de cambio (convergencia estavel, testei outros valores inic.)
%otimizador=[0.9 -1]
%objh = @(otim) ajuste(otim,set,IRFrex,pinv(Wrex),mod,'epsrstar')/1000-14;
%[param_est,wlf_est] = fmincon(objh, otimizador,[],[],[],[],[0,-100,1.01],[1,0,10],[],options);
%para o choque monetario (convergencia estavel, testei outros valores inic.)
%otimizador=[0.55 8]
%ajuste(otimizador,set,IRFr,pinv(Wr),mod,'epsmm')/1000-3.26;
%objh = @(otim) ajuste(otim,set,IRFr,pinv(Wr),mod,'epsmm')/1000-3.26;
%[param_est,wlf_est] = fmincon(objh, otimizador,[],[],[],[],[0,0],[1,100],[],options);
%

