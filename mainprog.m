clear all
close all

%%**********************************************************
% Resolve o modelo analiticamente
%**********************************************************
%parametros
[param0,set] = parameters_est;
%priors=prior(param0) %priors em formato structure (mudei para array para agilizar no loop do mcmc, entao essa funcao ficou meio inutil

%Define o modelo e calcula as matrizes do gensys analiticamente
model = modelsims(param0,set);
%dimensao do modelo
neq=length(model.f);
%Gerar uma funcao que retorna as matrizes do Gensys 
%essa linha criara um file ".m" com as matrizes simbolicas que
%serao usadas na solucao do modelo (nome do arquivo criado: model_prog)
%aumentar eficiencia: nao precisa computar as derivadas em toda iteracao
model_func_sims(model);

%Transforma em matriz os parametros fixos
set = struct2array(set);

%**********************************************************
%% Selecao das observaveis e matriz de selecao do KF
%**********************************************************
variaveis={'pidd' 'pif' 'dy' 'xl'  'pix'  'r'};

%selecinona as observaveis
[dados txt raw]=xlsread('data.xlsx','data_77','B1:I119');
for i=1:length(variaveis);
	idx(i)=find(strcmp(txt,variaveis{i}));
end
dados=dados(:,[idx]);

%Matriz de selecao
for i=1:length(variaveis)
    idx(i)=find(model.Y==variaveis{i});
end
H=zeros(length(idx),neq); 
for i=1:length(idx);
	H(i,idx(i))=1;
end
clear idx;

%% Kalman filter

warm=0;  %quantidade de dados da amostra que serao descartados para computar a verossimilhanca dos dados
shat=zeros(neq,1);  %chute inicial dos estados
sig=eye(neq);    %chute inicial da matriz de cov.

%*************************************************************************
%*********************MAXIMIZACAO DA VEROSSIMILHANCA**********************
%*************************************************************************

%%Calcula valores inicias para os estados com base nas priors
%%(Lembrando: os estados sao transformacoes dos parametros !) 
sig_mode_mcmc;

%media dos states
 load('mean_prior');
% 
% %funcao objetivo
 obj=@(STATE) -(kf_ll(dados',state2param(STATE),set,warm,H,sig,shat)+logpriorMCMC(STATE,state2param(STATE)))
% 
%hessiana no ponto inicial
H0=numhessian(obj,mean0);
save('H0','H0');
load('H0');

%minimizacao codigo do SIMS
[fh, x0_LL, gh, sigma0_LL, itct]=csminwel(obj,mean0,H0^(-1),[],1e-9,100000);    
%salva os resultados da maximizacao por ver

maxLL.fh=fh; maxLL.x0=x0_LL; maxLL.sigma0=sigma0_LL; maxLL.itct=itct;
save('max_LL','max_LL');
mode_LL=state2param(x0_LL);
save('mode_LL','mode_LL')

%% %***********************************************************************
%*************************************************************************
%*********************MARKOV CHAIN MONTE CARLO*****************************
%*************************************************************************
%*************************************************************************

%%Calcula valores inicias para os estados com base nas priors
%%(Lembrando: os estados sao transformacoes dos parametros !) 
sig_mode_mcmc;


%Load dos valoes iniciais (baseados nas priors)
load('sigma_prior') ;  %chute inicial para a matriz de varianca da cadeia de markov do MCMC
load('mode_prior')   ;%ponto inicial para a cadeia de markov do MCMC
load('mean_prior');

%%define a funcao que calcula a densidade dos dados
   ll_kf=@(PAR) kf_ll(dados',PAR,set,warm,H,sig,shat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parte 1: Estimacao preliminar da moda e covariancia da posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %DESATIVAR ESSE PASSO SE FOR USAR OS VALORES OBTIDOS NA MAXIMIZACAO
% 
 tun=1;
 ndraws=5000;
 nburn=0.4*ndraws;
 sigma0=diag(sigma_prior)
 x0=mean_prior;
 %profile on
     [x accept_rate mode_state]=mcmc(ndraws,x0,sigma0,tun,ll_kf);
 %profile viewer
 
 %Burna e calcula moda e covariancia
 x=x(nburn:end,:);
 sigma0=cov(x);
 x0=mode_state;
 
 for i=1:size(x,1);
 param(i,:)= state2param(x(i,:));
 end
 
 save('sigma0','sigma0');save('x0','x0');
 
 %salva em um arquivo diferente para comparar depois com as obitdas na LL
 sigma0_MC=sigma0; x0_MC=x0;
 save('sigma0_MC','sigma0');save('x0_MC','x0');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parte 2: Tunning para o acceptance ratio e matriz de variancia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Usa moda e variancia obtidas pelo processo de monte carlo

%     load('sigma0','sigma0');load('x0','x0');%load('sig','sig');load('shat','shat');
% 
% %Usa moda e varianca obtidas da maximizacao da verossimilhanca
% % Note que o csminwel MINIMIZA. Entao o inverso da hessiana no problema de
% % maximizacao teorico eh a propria matriz H que sai do csminwel.

     load('sigma0_LL'); load('x0_LL');
     sigma0=sigma0_LL; x0=x0_LL;

%inicializacao    
 tun=1;
 tun_=1;
 desired_rate=0.3;
 l_b=desired_rate-0.07; u_b=desired_rate+0.07;  %bounds para a aceitacao desejada
 max_attemp=100; %nro maximo de tentativas de tunnar 
 maxiter=1; %nro de repeticao do processo inteiro: vai achar "maxiter" tunning parameters otimos
 jump=2;   
 cunha=1000;
 
 for j=1:maxiter;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%2.1: Tunning
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ndraws=20000;
 for i=1:max_attemp
     [x accept_rate x0]=mcmc_busca_tun(ndraws,x0,sigma0,tun_,ll_kf,l_b,u_b);
     
     if (accept_rate)>u_b % se aceita demais
         if tun_>=tun  %testa se o tun anterior era menor que o atual. Se for, precisamos entao aumentar mais ainda o tun
             tun=tun_;
             tun_=tun_*jump;
         else            %caso o tun anterior ja fosse alto, pegamos a media
             tun_=(tun_+tun)/2;
         end
     elseif (accept_rate)<l_b %se aceita pouco
         if tun_<=tun  %testa se o tun atual e menor que o anterior. Se for e mesmo assim aceitamos demais, reduz mais ainda o tun
             tun=tun_;
             tun_=tun_/jump;
         else
             tun_=(tun_+tun)/2;
         end
     else %se caiu na regiao desejada,
         best_tun=tun_;
         tun=tun_;
         disp(strcat('busca ',num2str(j),' do parametro de tunning convergiu: acceptance rate= ',num2str(accept_rate)));
         break
     end
     
     %Verifica se o novo tun eh melhor que o anterior, guarda  e atualiza(se nao convergir alguma coisa ele pega)
     if abs(accept_rate-desired_rate)<cunha;
         cunha=abs(accept_rate-desired_rate);
         best_tun=tun;
     end
 i    
 end
 save('best_tun','best_tun'); save('x0','x0');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%2.2: Matriz de covariancia da proposed
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ndraws=20000;
 nburn=0.4*ndraws;
 
 [x accept_rate x0]=mcmc(ndraws,x0,sigma0,best_tun,ll_kf);
 
 %atualiza covariancia
 x=x(nburn:end,:);
 sigma0=cov(x);
 save('sigma0');
 %verifica se ao final do ultimo passo o tunning esta realmente no intervalo
 if j==maxiter
     if l_b<accept_rate & accept_rate<u_b
         disp(strcat('busca do parametro de tunning convergiu: acceptance ratio= ',num2str(accept_rate)));
     else
         prompt = strcat('Accept_rate=',num2str(accept_rate),' final fora do intervalo desejado. Deseja contiuar ? 1 para sim','Se sim, quantas iteracoes do processo inteiro ?');
         dlg_title = 'Input';
         num_lines = 1;
         def = {'1','3'};
         answer = inputdlg(prompt,dlg_title,num_lines,def);
         if eval(answer{1})==1;
             j=1; maxiter=eval(answer{2});
         end
     end
 end
 j
 end
 
 save('sigma0','sigma0');save('x0','x0');
 save('best_tun','best_tun');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parte 3: Posterior final
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('sigma0','sigma0');load('x0','x0');
load('best_tun');

disp('Calculando a posterior final');

%Posterior final: roda nblocks MCMC's com cadeias de markov comecando em
%pontos diferentes. Um choque normal e adicionado a moda calculada nos
%passos anteriores para alterar as condicoes iniciais.

ndraws=1000000;
nblocks=3;

X=zeros(ndraws,size(sigma0,1),nblocks);
ACCCEPT_RATE=zeros(nblocks,1);
SIGMA=zeros(size(sigma0,1),size(sigma0,2),nblocks);
mode_posteriors=zeros(length(x0),1,nblocks);
for j=1:nblocks
	load('x0','x0');
	x0=x0+mvnrnd(zeros(length(x0),1),best_tun*sigma0,1);
	%profile on
	[x accep_rate x0_]=mcmc(ndraws,x0,sigma0,best_tun,ll_kf);
	%profile viewer

	X(:,:,j)=x; SIGMA(:,:,j)=sigma0; ACCEPT_RATE(j)=accept_rate;
	mode_posteriors(:,:,j)=x0_;
	j
end

MCMC_LL.X=X; MCMC_LL.SIGMA=SIGMA; MCMC_LL.ACCEPT_RATE=ACCEPT_RATE; MCMC_LL.mode_posteriors=mode_posteriors;

save('!MCMC_LL','MCMC_LL');
