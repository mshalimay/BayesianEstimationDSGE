function [LL shat sig]=kf_ll(dados,PAR,set,warm,H,sig,shat)
%calcula a verossimilhanca do modelo atraves do filtro de kalman
%Inputs: 
%seja q = numero de variaveis no modelo q1=observaveis e T=tamanho da amostra
    %dados: matriz q1xT de dados
    %PAR: vetor de parametros a ser estimado
    %set: vetor de parametros fixos na estimacao
    %warm: percentual dos dados usados para aquecimento
    %H: matriz q1xq de selecao para as observaveis do conjunto de estados
    %shat: chute inicial para o vetor de estados (opcional)
    %sig: chute inicial para a matriz de covariancia dos estados (opcional)
%Outputs:
    % verossimilhanca dos dados dado o vetor de parametros
    
%Matrizes do Gensys    
[g0 g1 PSI PI]=model_prog(PAR,set);
C=zeros(size(g1,1),1);

%dimensoes
nvar=size(g1,1);
T=size(dados,2);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);

if eu(2)~=1
    %eu
    %PAR
    %error('Gensys sem solucao')
    LL=-Inf;  %nao pode usar -infinito por causa da minimizacao por verossimilhanca. Se nao for usa-la, troque para -Inf caso desjee
    shat=shat; sig=sig;
else

% switch nargin
%     case 5    
%         shat=zeros(nvar,1);  %chute inicial
%         sig=10*eye(nvar);    %
%     case 6
%         shat=zeros(nvar,1);  %chute inicial
% end

%% Aquecimento
if warm>0
    for i=1:warm;
        [shat,sig,lh,ydev]=kf(dados(:,i),H,shat,sig,G1,impact,C);
    end;
end;

LL=0;
%% Verossimilhanca
for i=(warm+1):T
    [shat,sig,lh,ydev]=kf(dados(:,i),H,shat,sig,G1,impact,C);
    LL=sum(lh)+LL;
end
end