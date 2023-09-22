function state2param = state2param(state)
%Essa funcao transforma os estados saindo da cadeia de markov no MCMC de
%volta nos parametros do modelo

%input: estados da cadeia de markov
%output: parametros do modelo associados

%% Transformacoes
%Os estados no MCMC sao draws de uma normal, mas param nao tera 
%necessariamente o mesmo suporte que o estado
%Solucao: gera os estados de uma transformacao dos parametros e retransforma o
%estado de volta pro parametro

%Seja theta o parametro do modelo e x o estado sorteado no MCMC
%Duas transformacoes sao feitas: t1 e t2.

%t1 garante que parametro sorteado no MCMC sao positivos: theta=exp(x)
%t2 garante que seja entre 0 e 1: theta=exp(x)/(1+exp(x))

state2param=zeros(length(state),1)+NaN;
t1=@(x) exp(x);
t2=@(x) exp(x)./(1+exp(x));
%cont=@(x) sum(isnan(x))
%% Calculo dos parametros
state2param(1)=t2(state(1)); %rigidez domestica
state2param(2)=t2(state(2)); %rigidez importados
state2param(3)=t2(state(3)); %rigidez exportados

state2param(4)=t2(state(4)); %indexacao domestica
state2param(5)=t2(state(5)); %indexacao importados
state2param(6)=t2(state(6)); %indexacao exportados

state2param(7)=state(7);     %elast. subst intertemp.
state2param(8)=state(8);     %elast. subst. importados x domesticos

%parametros da regra de taylor
state2param(9)=state(9);    %politica monetaria
state2param(10)=state(10);  %parametro do cresc. da iflacao na TR
state2param(11)=state(11);  %parametro do cambio na TR
state2param(12)=state(12); %parametro do cresc. pib na TR

%parametros autoregressivos dos processos exogenos 
state2param(13:17)=t2(state(13:17)); 

%desvio padrao dos choques exogenos
%Obs.: para a ivngamma temos beta=(alpha-1)*media, 
%onde:alpha eh o shape (degree freedom) e beta eh o scale
%priors de Smets Wouters 2006 e malind Adolfsson
state2param(18:end)=t1(state(18:end));

end