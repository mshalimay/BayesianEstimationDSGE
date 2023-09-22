function [mod] = model(param,set)

%parametros simbolicos
param_list = fieldnames(param);
syms(param_list{1:end});
for j = 1:length(param_list)
    eval(['PARAM(j) = ',param_list{j} ,';']);
end
PARAM

set_list = fieldnames(set);
syms(set_list{1:end});
for j = 1:length(set_list)
    eval(['SET(j) = ',set_list{j} ,';']);
end
%Variaveis simbolicas 
%Descr.: x_f eh a variavel no periodo seguinte, x_l_j eh a o lag j da variavel
%epsx sao os choques , etax sao os erros expectacionais
%expx eh o valor esperado em t de x
%expx_f	 eh o forward do valor esperado em t de x. So precisa pra manter as matrizes do Gensys quadradas (todas as derivadas vao ser zero)
syms pdstar r de c a y pidd pif pix pi veta rer xl cstar rstar choquem dy dpi tt ttx ttf choquepif choquepix choquec dc ...
    pidd_l pif_l pix_l pi_l r_l y_l c_l tt_l ttf_l ttx_l rer_l pdstar_l;

syms pdstar_f r_f de_f c_f a_f y_f pidd_f pif_f pix_f pi_f veta_f rer_f xl_f cstar_f rstar_f choquem_f dy_f dpi_f tt_f ttx_f ttf_f choquepif_f choquepix_f choquec_f dc_f ...
    pidd_l_f pif_l_f pix_l_f pi_l_f r_l_f y_l_f c_l_f tt_l_f ttf_l_f ttx_l_f rer_l_f pdstar_l_f;               

syms expde exppi expc exppidd exppif exppix;
syms expde_f exppi_f expc_f exppidd_f exppif_f exppix_f

syms etade etapi etac etapidd etapif etapix;
syms epsvet epscstar epsrstar epsmm epspif epspix epsc;


% Vetores de variaveis, choques expectactionais e choques exogenos
Y = [pdstar r de c a y pidd pif pix pi veta rer xl cstar rstar choquem dy dpi tt ttx ttf choquepif choquepix choquec dc ...
    pidd_l pif_l pix_l r_l y_l c_l pi_l tt_l ttf_l ttx_l rer_l pdstar_l expde exppi expc exppidd exppif exppix] ; %variaveis em t

Yf = [pdstar_f r_f de_f c_f a_f y_f pidd_f pif_f pix_f pi_f veta_f rer_f xl_f cstar_f rstar_f choquem_f dy_f dpi_f tt_f ttx_f ttf_f choquepif_f choquepix_f choquec_f dc_f ...
    pidd_l_f pif_l_f pix_l_f r_l_f y_l_f c_l_f pi_l_f tt_l_f ttf_l_f ttx_l_f rer_l_f pdstar_l_f expde_f exppi_f expc_f exppidd_f exppif_f exppix_f] ;  %t+1

et = [etade etapi etac etapidd etapif etapix]; %choques expectacionais

eps = [epsvet epscstar epsrstar epsmm epspif epspix epsc]; %choques exogenos

%Eqs. do modelo LINEARIZADAS
syms f
f(length(f)+1)=r-rstar_f-expde;  %UIP

f(length(f)+1)=-c-1/sigm*(r-exppi)+expc+(1-v)/v*phibar*a_f+choquec_f;

f(length(f)+1)=-pi+(1-alphac)*pidd+alphac*pif;

f(length(f)+1)=-pidd+bet/(1+gammad*bet)*exppidd+gammad/(1+gammad*bet)*pidd_l-(1-bet*xid)*(1-xid)/(xid*(1+gammad*bet))*...
    ((1-(1-alphay)*(1-alphac))*(tt)-(1-alphay)*(omeg*y+sigm*c)+(1-alphay)*(1+omeg)*veta_f);

f(length(f)+1)=-pif+bet/(1+gammaf*bet)*exppif+gammaf/(1+gammaf*bet)*pif_l-(1-bet*xif)*(1-xif)/(xif*(1+gammaf*bet))*(ttf-choquepif_f);

f(length(f)+1)=-pix+bet/(1+gammax*bet)*exppix+gammax/(1+gammax*bet)*pix_l-(1-bet*xix)*(1-xix)/(xix*(1+gammax*bet))*(ttx-choquepix_f);

f(length(f)+1)=-y+(1-alphay)*(1-alphac)*c+(1-(1-alphay)*(1-alphac))*(cstar_f-eta*pdstar)-eta*alphac*(1-alphay)*(1-alphac)*(tt);

f(length(f)+1)=-a_f+rrbar*(a+y-c+(alphac+alphay/(1-alphay))*(-ttx));

f(length(f)+1)=-r+rhor*r_l+(1-rhor)*(phi*pi+phiy*y+phirer*rer)+rdpi*dpi+rdy*dy+choquem_f; %Regra de taylor

%processos exogenos
f(length(f)+1)=-veta_f+rhoveta*veta+sdvet*epsvet;
f(length(f)+1)=-cstar_f+rhocstar*cstar+sdcstar*epscstar;
f(length(f)+1)=-rstar_f+rhorstar*rstar+sdrstar*epsrstar;
f(length(f)+1)=-choquem_f+rhom*choquem+sdm*epsmm;
f(length(f)+1)=-choquepif_f+rhopif*choquepif+sdpif*epspif;
f(length(f)+1)=-choquepix_f+rhopix*choquepix+sdpix*epspix;
f(length(f)+1)=-choquec_f+rhoc*choquec+sdc*epsc;
%law of movement
f(length(f)+1)=-rer+rer_l+pi-de;
f(length(f)+1)=-tt+tt_l+pidd-pif;
f(length(f)+1)=-ttf+ttf_l+pif-de;
f(length(f)+1)=-ttx+ttx_l+pix-pidd+de;
f(length(f)+1)=-pdstar+pdstar_l+pix;

%ident
f(length(f)+1)=-xl+y-(1-alphac)*(1-alphay)*c;
f(length(f)+1)=-dpi+pi-pi_l;
f(length(f)+1)=-dy+y-y_l;
f(length(f)+1)=-dc+c-c_l;

%Equacoes de lag
f(length(f)+1)=-pidd_l_f+pidd;
f(length(f)+1)=-pif_l_f+pif;
f(length(f)+1)=-pix_l_f+pix;
f(length(f)+1)=-pi_l_f+pi;
f(length(f)+1)=-r_l_f+r;
f(length(f)+1)=-y_l_f+y;
f(length(f)+1)=-c_l_f+c;
f(length(f)+1)=-tt_l_f+tt;
f(length(f)+1)=-ttx_l_f+ttx;
f(length(f)+1)=-ttf_l_f+ttf;
f(length(f)+1)=-rer_l_f+rer;
f(length(f)+1)=-pdstar_l_f+pdstar;


%eqs expectacionais
f(length(f)+1)=-exppi+pi_f+etapi;
f(length(f)+1)=-exppidd+pidd_f+etapidd;
f(length(f)+1)=-exppif+pif_f+etapif;
f(length(f)+1)=-exppix+pix_f+etapix;
f(length(f)+1)=-expc+c_f+etac;
f(length(f)+1)=-expde+de_f+etade;

%remove colunas zeradas da matriz de eqs.
f(1)=0;
f(all(isAlways(f==0) ,2) ,:) = [];
f(:,all( isAlways(f==0) ,1)) = [];

mod.Y = Y;
mod.Yf = Yf;
mod.eps = eps;
mod.et = et;
mod.PARAM = PARAM;
mod.param = param;
mod.SET = SET;
mod.set = set;

mod.f = f;
mod = anal_deriv_sims(mod);


