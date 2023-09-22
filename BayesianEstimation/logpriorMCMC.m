function [prior] = logpriorMCMC(state,param)
%Essa funcao avaliara a prior em cada parametro saindo dos draws do MCMC

%input: Vetor de estados da cadeia de markov do MCMC e parametros associados (ja transformados)
%output: vetor de priors dos estados

priors=zeros(length(state),1);

%% Transformacoes de variaveis no MCMC.
%Os estados no MCMC sao draws de uma normal, mas param nao tera necessariamente o mesmo suporte que o estado
%Solucao: gera os estados de uma variavel transformada e retransforma o
%estado de volta pro parametro
%funcao state2param transforma os estados da cadeia de markov de volta nos
%parametros do modeloo

%Alem de fazer o mapping de volta para os parametros, eh preciso ajustar a 
%proposal pelo jacobiano da transformacao feita

%Ao final do MCMC a distribuicao recuperada sera a distribuicao dos
%parametros TRANSFORMADOS (ou seja, dos estados do MCMC). eh preciso
%transformar de volta para obter a densidade dos parametros do modelo

%Seja theta o parametro como no modelo e x o estado sorteado no MCMC
%Duas transformacoes sao feitas: t1 e t2.

%t1 garante que parametro sorteado no MCMC sao positivos: theta=exp(x)
%t2 garante que seja entre 0 e 1: theta=exp(x)/(1+exp(x))
%j1 e j2 sao os jacobianos das transformacoes acima no formato log.
%j1=exp(x) --> log(j1)=x;  
%j2=exp(x)/(1+exp(x))^2 --> log(j2)=x-2*log(1+exp(x)))

%A posterior final eh dada por: p(x|Y)~p(Y|t(x))*prior(t(x))*j

%j1=@(x) x;
j2=@(x) x-2*log((1+exp(x)));

%% Calculo das priors
%rigidez
priors(1)=log(beta_dens(param(1),0.675,0.05))+j2(state(1)); %beta media 0.675 sd 0.05
priors(2)=log(beta_dens(param(2),0.5,0.1))+j2(state(2)); %beta media 0.5 sd 0.1
priors(3)=log(beta_dens(param(3),0.5,0.1))+j2(state(3)); %beta media 0.5 sd 0.1
%indexacao
priors(4)=log(beta_dens(param(4),0.5,0.15))+j2(state(4)); %beta media 0.5 e sd 0.15
priors(5)=log(beta_dens(param(5),0.5,0.15))+j2(state(5));%beta media 0.5 e sd 0.15
priors(6)=log(beta_dens(param(6),0.5,0.15))+j2(state(6));%beta media 0.5 e sd 0.15

priors(7)=log(normpdf(param(7),1.5,0.375));
priors(8)=log(normpdf(param(8),1.5,0.5));

%parametros da regra de taylor
priors(9)=log(normpdf(param(9),1.7,0.1));%*(param(9)>=1.1 & param(9)<100);
priors(10)=log(normpdf(param(10),0.3,0.1));
priors(11)=log(normpdf(param(11),0,0.05));
priors(12)=log(normpdf(param(12),0.0625,0.05));

%parametros autoregressivos dos processos exogenos 
priors(13)=log(beta_dens(param(13),0.9,0.1))+j2(state(13)); %beta media 0.9, sd 0.1
priors(14)=log(beta_dens(param(14),0.9,0.1))+j2(state(14)); %beta media 0.9, sd 0.1
priors(15)=log(beta_dens(param(15),0.9,0.1))+j2(state(15)); %beta media 0.9, sd 0.1
priors(16)=log(beta_dens(param(16),0.9,0.1))+j2(state(16)); %beta media 0.9, sd 0.1
priors(17)=log(beta_dens(param(17),0.9,0.1))+j2(state(17)); %beta media 0.9, sd 0.1


%desvio padrao dos choques exogenos
%Obs.: para a ivngamma temos beta=(alpha-1)*media, 
%onde:alpha eh o shape (degree freedom) e beta eh o scale
%priors de Smests Wouters 2006
priors(18)=log(invgamma_dens(param(18),0.15,Inf))+state(18);
priors(19)=log(invgamma_dens(param(19),0.05,Inf))+state(19);
priors(20)=log(invgamma_dens(param(20),0.7,Inf))+state(20);
priors(21)=log(invgamma_dens(param(21),0.3,Inf))+state(21);
priors(22)=log(invgamma_dens(param(22),0.3,Inf))+state(22);
priors(23)=log(invgamma_dens(param(23),0.2,Inf))+state(23);
priors(24)=log(invgamma_dens(param(24),0.2,Inf))+state(24);

prior=sum(priors);
end