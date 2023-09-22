
%Esse codigo gera uma estimativa inicial para as variancias e moda
%dos estados do MCMC pq nao podemos usar os da prior diretamente pois
%os estados da cadeia de markov sao transf. dos parametros originais

%A relacao entre os parametros do modelo e os estados da cadeia de markov
%eh dada por:
%theta e o parametro; x e o state
%Para theta>0: theta=exp(x)
%Para 0<theta<1: theta=exp(x)/(1+exp(x))

%Portanto, a transformacao inversa para obter os estados eh dada por:
% para theta>0: x=log(theta)
% para 0<theta<1: x = log (theta/(1-theta))

%% Calculo das variancias e modas
%Obs.: criei uma funcao que calcula a moda pq a do matlab eh uma merda
sigma_prior=zeros(24,1);
mode_prior=zeros(24,1);
mean_prior=zeros(24,1);

%rigidez domestica
x=betarnd(58.5563,28.1937,1000000,1); %beta media 0.675 sd 0.05
x=log(x./(1-x)); 
sigma_prior(1)=var(x);
mode_prior(1)=moda(x);
mean_prior(1)=mean(x);

%rigidez importados e exportados
x=betarnd(12,12,1000000,1); %beta media 0.5 sd 0.1
x=log(x./(1-x)); 
sigma_prior(2:3)=var(x);
mode_prior(2:3)=moda(x);
mean_prior(2:3)=mean(x);

%indexacoes
x=betarnd(5.05556,5.05556,1000000,1); %beta media 0.5 sd 0.15
x=log(x./(1-x)); 
sigma_prior(4:6)=var(x);
mode_prior(4:6)=moda(x);
mode_prior(4:6)=mean(x);

%estado e parametro sao normais(x,sig^2)
sigma_prior(7)=0.375^2;  %elast. subst. intertemporal
mode_prior(7)=1.5;
mean_prior(7)=1.5;

sigma_prior(8)=0.5^2;    %elast. subst. importados exportados
mode_prior(8)=1.5;
mean_prior(8)=1.5;

%parametros da regra de taylor
sigma_prior(9)=0.1^2;   %inflation response
mode_prior(9)=1.7;
mean_prior(9)=1.7;

sigma_prior(10)=0.1^2; %inflation growth response
mode_prior(10)=0.3;
mean_prior(10)=0.3;

sigma_prior(11)=0.05^2; %rer response
mode_prior(11)=0;
mean_prior(11)=0;

sigma_prior(12)=0.05^2;   %gdp growth response
mode_prior(12)=0.0625;
mean_prior(12)=0.0625;

%parametros autoregressivos dos processos exogenos 
x=betarnd(7.2,0.8,1000000,1); %beta media 0.9 sd 0.1
x=log(x./(1-x)); 
sigma_prior(13:17)=var(x);
mode_prior(13:17)=moda(x);
mean_prior(13:17)=mean(x);

%desvio padrao dos choques exogenos
%epsm
x=1./gamrnd(2,1/((2-1)*0.15),10000000,1); % invgamma media=0.15 e shape=degrees freedom=2
x=log(x);
sigma_prior(18)=var(x);
mode_prior(18)=moda(x);
mean_prior(18)=mean(x);

%epsrstar
x=1./gamrnd(2,1/((2-1)*0.05),10000000,1); % invgamma media=0.05 e shape=degrees freedom=2
x=log(x);
sigma_prior(19)=var(x);
mode_prior(19)=moda(x);
mean_prior(19)=mean(x);

%epsveta
x=1./gamrnd(2,1/((2-1)*0.7),10000000,1); % invgamma media=0.7 e shape=degrees freedom=2
x=log(x);
sigma_prior(20)=var(x);
mode_prior(20)=moda(x);
mean_prior(20)=mean(x);

%epspif e epspix
x=1./gamrnd(2,1/((2-1)*0.3),10000000,1); % invgamma media=0.3 e shape=degrees freedom=2
x=log(x);
sigma_prior(21:22)=var(x);
mode_prior(21:22)=moda(x);
mean_prior(21:22)=mean(x);

%epscstar e epsc
x=1./gamrnd(2,1/((2-1)*0.2),10000000,1); % invgamma media=0.2 e shape=degrees freedom=2
x=log(x);
sigma_prior(23:24)=var(x);
mode_prior(23:24)=moda(x);
mean_prior(23:24)=mean(x);


save('sigma_prior','sigma_prior')
save('mode_prior','mode_prior')
save('mean_prior','mean_prior')
