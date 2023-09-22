clear
close all

%%
%**********************************************************
% Resolve o modelo analiticamente - sem optimal policy
%**********************************************************

%parametros
[param0,setted] = parameters_sim;
%Define o modelo e calcula as matrizes do gensys analiticamente
mod = modelsims(param0,setted);
%dimensao do modelo
neq=length(mod.f);
%Gerar uma funcao que retorna as matrizes do Gensys 
%essa linha criara um file ".m" com as matrizes simbolicas que
%serao usadas na solucao do modelo (nome do arquivo criado: model_prog)
%aumentar eficiencia: nao precisa computar as derivadas em toda iteracao
model_func_sims(mod);
var={'y' 'c' 'xl' 'tt' 'rr' 'rer' 'pidd' 'pif' 'w' 'a'};
eps={'epsvet','epscstar','epsrstar'};
idxv=zeros(1,size(var,2));
idxeps=zeros(1,3)

for i=1:size(var,2);
   idxv(i)=find(mod.Y==var{i});  
end
for i=1:3;
   idxeps(i)=find(mod.eps==eps{i});  
end

setted = struct2array(setted);

%%
%*************************************************************************
%Calculo da funcao de resposta ao impulso - TR na inflacao domestica
%*************************************************************************
%Valores iniciais dos parametros
param = [0.85,0.85,0,0];
%Matrizes do Gensys com valores iniciais
[g0 g1 PSI PI]=model_prog(param,setted);
C=zeros(neq,1);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);

T=100;
IRv=zeros(T,neq);IRrstar=zeros(T,neq);IRcstar=zeros(T,neq);IRm=zeros(T,neq);
%Ajustando os choques para coincidir com as dinamicas do paper
tunning=ones(4,1);
tunning=zeros(size(impact,1),1);
tunning(find(mod.eps=='epsvet'))=1.25;
tunning(find(mod.eps=='epsrstar'))=-0.45; %para bater com os ultimos graf.
tunning(find(mod.eps=='epscstar'))=2.4;
tunning(find(mod.eps=='epsmm'))=1; 
%Calcula a impulso resposta
[IR]=ir_gensys(G1,impact,T,tunning) ;%tunning opcional para ajustar o tamanho do choque na IRF. Se vazio, tunning=1
%transforma Impulsos de inflacao e juros em anuais (como no paper)
IR(:,find(mod.Y=='pidd'),:)=4*IR(:,find(mod.Y=='pidd'),:);
IR(:,find(mod.Y=='pif'),:)=4*IR(:,find(mod.Y=='pif'),:);
IR(:,find(mod.Y=='pi'),:)=4*IR(:,find(mod.Y=='pi'),:);
IR(:,find(mod.Y=='rr'),:)=4*IR(:,find(mod.Y=='rr'),:);
% %coleta impulsos de cada choque em vetores especificos
% IRrstar=IR(:,:,find(mod.eps=='epsrstar'));
% IRv=IR(:,:,find(mod.eps=='epsvet'));
% IRcstar=IR(:,:,find(mod.eps=='epscstar'));
% IRm=IR(:,:,find(mod.eps=='epsmm'));

IRtr=zeros(T,size(idxv,2),size(idxeps,2))

for j=1:size(idxeps,2)
    for i=1:size(idxv,2)
        IRtr(:,i,j)=IR(:,idxv(i),idxeps(j))
    end
end
%%
%*************************************************************************
%Calculo da funcao de resposta ao impulso - Precos flexiveis
%*************************************************************************
%Valores iniciais dos parametros
param = [0.0001,0.0001,0,0];
%Matrizes do Gensys com valores iniciais
[g0 g1 PSI PI]=model_prog(param,setted);
C=zeros(neq,1);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);

T=100;
IRv=zeros(T,neq);IRrstar=zeros(T,neq);IRcstar=zeros(T,neq);IRm=zeros(T,neq);
%Ajustando os choques para coincidir com as dinamicas do paper
tunning=ones(4,1);
tunning=zeros(size(impact,1),1);
tunning(find(mod.eps=='epsvet'))=1.25*10^8;
tunning(find(mod.eps=='epsrstar'))=-0.45; %para bater com os ultimos graf.
tunning(find(mod.eps=='epscstar'))=2.4;
tunning(find(mod.eps=='epsmm'))=1; 
%Calcula a impulso resposta
[IR]=ir_gensys(G1,impact,T,tunning) ;%tunning opcional para ajustar o tamanho do choque na IRF. Se vazio, tunning=1
%transforma Impulsos de inflacao e juros em anuais (como no paper)
IR(:,find(mod.Y=='pidd'),:)=4*IR(:,find(mod.Y=='pidd'),:);
IR(:,find(mod.Y=='pif'),:)=4*IR(:,find(mod.Y=='pif'),:);
IR(:,find(mod.Y=='pi'),:)=4*IR(:,find(mod.Y=='pi'),:);
IR(:,find(mod.Y=='rr'),:)=4*IR(:,find(mod.Y=='rr'),:);
%coleta impulsos de cada choque em vetores especificos
% IRrstar=IR(:,:,find(mod.eps=='epsrstar'));
% IRv=IR(:,:,find(mod.eps=='epsvet'));
% IRcstar=IR(:,:,find(mod.eps=='epscstar'));
% IRm=IR(:,:,find(mod.eps=='epsmm'));

IRflex=zeros(T,size(idxv,2),size(idxeps,2))

for j=1:size(idxeps,2)
    for i=1:size(idxv,2)
        IRflex(:,i,j)=IR(:,idxv(i),idxeps(j))
    end
end

%%
%%
%**********************************************************
% Resolve o modelo analiticamente para Optimal Policy
%**********************************************************

%parametros
[param,setted] = parameters_sim;
%Define o modelo e calcula as matrizes do gensys analiticamente
mod = model_sims_opt(param,setted);
%dimensao do modelo
neq=length(mod.f);
%Gerar uma funcao que retorna as matrizes do Gensys 
%essa linha criara um file ".m" com as matrizes simbolicas que
%serao usadas na solucao do modelo (nome do arquivo criado: model_prog)
%aumentar eficiencia: nao precisa computar as derivadas em toda iteracao
model_func_sims(mod);
var={'y' 'c' 'xl' 'tt' 'rr' 'rer' 'pidd' 'pif' 'w' 'a'};
eps={'epsvet','epscstar','epsrstar'};
idxv=zeros(1,size(var,2));
idxeps=zeros(1,3)

for i=1:size(var,2);
   idxv(i)=find(mod.Y==var{i});  
end
for i=1:3;
   idxeps(i)=find(mod.eps==eps{i});  
end

setted = struct2array(setted);


%*************************************************************************
%Calculo da funcao de resposta ao impulso - Optimal policy
%*************************************************************************
%Valores iniciais dos parametros
param = [0.85,0.85,0,0];
%Matrizes do Gensys com valores iniciais
[g0 g1 PSI PI]=model_prog(param,setted);
C=zeros(neq,1);

%Gensys
[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(-g0,g1,C,PSI,PI);

T=100;
IRv=zeros(T,neq);IRrstar=zeros(T,neq);IRcstar=zeros(T,neq);IRm=zeros(T,neq);
%Ajustando os choques para coincidir com as dinamicas do paper
tunning=ones(4,1);
tunning=zeros(size(impact,1),1);
tunning(find(mod.eps=='epsvet'))=1.25;
tunning(find(mod.eps=='epsrstar'))=-0.45; %para bater com os ultimos graf.
tunning(find(mod.eps=='epscstar'))=2.4;
tunning(find(mod.eps=='epsmm'))=1; 
%Calcula a impulso resposta
[IR]=ir_gensys(G1,impact,T,tunning) ;%tunning opcional para ajustar o tamanho do choque na IRF. Se vazio, tunning=1
%transforma Impulsos de inflacao e juros em anuais (como no paper)
IR(:,find(mod.Y=='pidd'),:)=4*IR(:,find(mod.Y=='pidd'),:);
IR(:,find(mod.Y=='pif'),:)=4*IR(:,find(mod.Y=='pif'),:);
IR(:,find(mod.Y=='pi'),:)=4*IR(:,find(mod.Y=='pi'),:);
IR(:,find(mod.Y=='rr'),:)=4*IR(:,find(mod.Y=='rr'),:);
% %coleta impulsos de cada choque em vetores especificos
% IRrstar=IR(:,:,find(mod.eps=='epsrstar'));
% IRv=IR(:,:,find(mod.eps=='epsvet'));
% IRcstar=IR(:,:,find(mod.eps=='epscstar'));
% IRm=IR(:,:,find(mod.eps=='epsmm'));

IRopt=zeros(T,size(idxv,2),size(idxeps,2))

for j=1:size(idxeps,2)
    for i=1:size(idxv,2)
        IRopt(:,i,j)=IR(:,idxv(i),idxeps(j))
    end
end
%%
%*************************************************************************
%Plot dos graficos
%*************************************************************************
titulos={'output','consumption','net export','terms of trade','real interest rate', 'real exchange rate', 'inflation', 'import inflation', 'real wage', 'net foreign assets'}
t=20;
zero=zeros(1,t)
for j=1:size(idxeps,2)
   fig=figure
   color=get(fig,'Color')
   for i=1:size(idxv,2)
        subplot(5,2,i)
        if strcmp(titulos{i},'inflation')==1 | strcmp(titulos{i},'import inflation')==1
            h=plot(1:t,IRtr(1:t,i,j),'red',1:t,IRopt(1:t,i,j),'blue');
            set(h(1),'LineStyle','--','LineWidth',1.5)
        else
            h=plot(1:t,IRtr(1:t,i,j),'red',1:t,IRflex(1:t,i,j),'green',1:t,IRopt(1:t,i,j),'blue');
            set(h(2),'LineWidth',1)
            set(h(1),'LineStyle','--','LineWidth',1.5)
        end
        title(titulos(i))
        box off
        hold on
        hline(0,'black')
        hold off
        if i==1 
            legend('Taylor Rule','Flex','Optimal','Position',[0.5,0.5,420,25],'Orientation','Horizontal')
        end
   end
saveas(fig,strcat(eps{j},'.png')) 
end
