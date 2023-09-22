function [x accept_rate mode_state]= mcmc(ndraws,x0,sigma,tun,ll_kf)

x=zeros(ndraws,length(x0));
%x(1,:)=mvnrnd(x0,tun*sigma);

%condicoes iniciais
x(1,:)=x0;
par=state2param(x(1,:));
post=ll_kf(par)+logpriorMCMC(x(1,:),par);
%choques da cadeia de markov
e=mvnrnd(zeros(length(par),1),tun*sigma,ndraws);
%probabilidades de aceitacao
u=log(rand(ndraws,1));

%inicializacao
reject=0;
accept=0;
max_post=post;
mode_state=x0;

%parametros para estimacao iterativa da variancia da posterior (como no
%mode_compute=6 do dynare)
%mi=x0';
%var=sigma;
%t=0;

for i=2:(ndraws)
   %t=t+1;
   x_=x(i-1,:)+e(i,:);
   par=state2param(x_);
   
   [T,llkf]=evalc('ll_kf(par)');
   post_=llkf+logpriorMCMC(x_,par);

   %teste 1: aceita ou rejeita o novo ponto
   if u(i)<post_-post  
       x(i,:)=x_;  
       %Teste 2: muda ou nao a moda
       if post_>max_post;
            max_post=post_;
            mode_state=x_;
       end
       post=post_;
       accept=accept+1;
   else %rejeitou
        x(i,:)=x(i-1,:);
       reject=reject+1;
   end
   %mi_=mi+1/t*(x(i,:)-x(i-1,:));
  
   %var_= var + mi'*mi - mi_'*mi_ + 1/t*(x(i,:)'*x(i,:)-var-mi'*mi);
   
   %mi=mi_;var=var_;
   if mod(i,500)==0
        i
        accept_rate=accept/(reject+accept)
   end
end
accept_rate=accept/(reject+accept);




    







