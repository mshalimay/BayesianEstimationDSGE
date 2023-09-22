var pdstar r de c a y pidd pif pix pi veta rer xl cstar rstar choquem dy dpi tt ttx ttf choquepif choquepix choquec dc;
varexo epsveta epscstar epsrstar epsm epspif epspix epsc;
parameters beta alphay alphac v rrbar phibar xid xif xix gammad gammaf gammax sigma omega eta rhor phi rdpi phirer phiy rdy rhoveta rhorstar rhocstar rhom
 rhopif rhopix rhoc;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %OBS.: identification is highly sensible to v
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%fixed parameters
beta=0.99;
alphay=0;
alphac=0.3;
v=.85;
rrbar=1/beta;
phibar=1-v*beta;
%rhom=0;
%sigma=1;
omega=2;
%eta=5;
rhor=0.8;
phiy=0.125;
rhoc=0.9;
rhopif=0.9;

model(linear);
r=rstar+de(+1);

c=-1/sigma*(r-pi(+1))+c(+1)+(1-v)/v*phibar*a+choquec;

pi=(1-alphac)*pidd+alphac*pif;

pidd=beta/(1+gammad*beta)*pidd(+1)+gammad/(1+gammad*beta)*pidd(-1)-(1-beta*xid)*(1-xid)/(xid*(1+gammad*beta))*
((1-(1-alphay)*(1-alphac))*(tt)-(1-alphay)*(omega*y+sigma*c)+(1-alphay)*(1+omega)*veta);

pif=beta/(1+gammaf*beta)*pif(+1)+gammaf/(1+gammaf*beta)*pif(-1)-(1-beta*xif)*(1-xif)/(xif*(1+gammaf*beta))*(ttf-choquepif);

	pix=beta/(1+gammax*beta)*pix(+1)+gammax/(1+gammax*beta)*pix(-1)-(1-beta*xix)*(1-xix)/(xix*(1+gammax*beta))*(ttx-choquepix);

y=(1-alphay)*(1-alphac)*c+(1-(1-alphay)*(1-alphac))*(cstar-eta*pdstar)-eta*alphac*(1-alphay)*(1-alphac)*(tt);

a=rrbar*(a(-1)+y-c+(alphac+alphay/(1-alphay))*(-ttx));

r=rhor*r(-1)+(1-rhor)*(phi*pi+phiy*y+phirer*rer)+rd	pi*dpi+rdy*dy+choquem; %Taylor rule 

%exogenous process
veta=rhoveta*veta(-1)+epsveta;
cstar=rhocstar*cstar(-1)+epscstar;
rstar=rhorstar*rstar(-1)+epsrstar;
choquem=rhom*choquem(-1)+epsm;
choquepif=rhopif*choquepif(-1)+epspif;
choquepix=rhopix*choquepix(-1)+epspix;
choquec=rhoc*choquec(-1)+epsc;

%law of movement
tt=tt(-1)+pidd-pif;
ttf=ttf(-1)+pif-de;
ttx=ttx(-1)+pix-pidd+de;
rer=rer(-1)+pi+de;
pdstar=pdstar(-1)+pix;
%identities
xl=y-(1-alphac)*(1-alphay)*c;
dpi=pi-pi(-1);
dy=y-y(-1);
dc=c-c(-1);
end;

steady_state_model;
 pdstar	  =0;
 r        =0;
 de       =0;
 c        =0;
 a        =0;
 y        =0;
 pidd      =0;
 pif      =0;
 pix      =0;
 pi       =0;
 veta     =0;
 rer      =0;
 xl       =0;
 cstar    =0;
 rstar    =0;
 choquem  =0;
 dy       =0;
 dpi      =0;
 tt       =0;
 ttx      =0;
 ttf      =0;
 choquepif =0; 
 choquepix =0; 
 choquec   =0;
 dc=0;
 end;

estimated_params;
xid, beta_pdf, 0.675, 0.05; %domestic rigidity
xif, beta_pdf, 0.5, 0.1;    %import rigidity
xix, beta_pdf, 0.5, 0.1;    %export rigidity

gammad, beta_pdf, 0.5, 0.15;	%domestic indexation
gammaf, beta_pdf, 0.5,0.15;		%import
gammax, beta_pdf, 0.5,0.15;		%export

sigma, normal_pdf, 1.5, 0.375;   %elast. intertemp. subst.
%omega, normal_pdf, 2, 0.75;	 	%elast. frisch
eta, normal_pdf, 1.5, 0.5;		%elast. subst. import x domestic goods

%rhor, beta_pdf, 0.8, 0.05; 		%interest rate smoothing taylor rule
phi, normal_pdf, 1.7, 0.1;		%inflation resp	onse TR
rdpi, normal_pdf, 0.3, 0.1;  	%inflation growth response TR
phirer, normal_pdf, 0, 0.05; 	%real exchange rate response TR
%phiy, normal_pdf, 0.125, 0.05; 	%gdp response TR
rdy, normal_pdf, 0.0625, 0.05;  %gdp growth response TR

rhoveta, beta_pdf, 0.9, 0.1;	  %productivity shock persistence
rhorstar, beta_pdf, 0.9, 0.1; %exchange rate shock persistence
rhocstar, beta_pdf, 0.9, 0.1;	  %external demand
rhom, beta_pdf, 0.9, 0.1; 	%monetary shock
%rhopif,beta_pdf,0.9,0.1;		%mark-up import sector
rhopix,beta_pdf,0.9,0.1;		%mark-up export sector
%rhoc,beta_pdf,0.9,0.1;			%preference shock

%inv gamma mean x and 2 degree freedom (shape) (or sd = infinity)
stderr epsm, inv_gamma_pdf, 0.45/3, inf; 		 %moetary shock 
stderr epsrstar, inv_gamma_pdf, 0.15/3, inf;   %exchange rate shock 
stderr epsveta, inv_gamma_pdf, 2.1/3, inf;     %product. shock
stderr epspif, inv_gamma_pdf, 0.9/3, inf;     %mark-up import
stderr epspix, inv_gamma_pdf, 0.9/3, inf;     %mark-up export
stderr epscstar, inv_gamma_pdf, 0.6/3, inf;   %foreign demand shock
stderr 	epsc, inv_gamma_pdf, 0.6/3, inf;       %preference shock
end;

varobs pidd pif dy xl pix	r;

 
%estimated_params_init(USE_CALIBRATION);	
// xid, 0.675;
// xif, 0.5;
// xix, 0.5;
// gammad,0.5;
// gammaf, 0.5;		%indexacao importados
// gammax, 0.5;		%indexacao exportados

// sigma, 1.5;   %elast. substituicao intertemp.
// omega, 2;	 	%elast. frisch
// eta,1.5;		%elast. subst. importados domesticos

// rhor,0.8; 		%parametro autoregressivo da TR
// phi,1.7;		%inflacao na TR
// rdpi,0.3;  	%parametro do cresc. da iflacao na TR
// phirer,0; 	%parametro do cambio na TR
// phiy,0.125; 	%parametro do pib na TR
// rdy,0.0625;  %parametro do cresc. pib na TR

// rhoveta,0.85;	%persistencia choque produtiv.
// rhorstar,0.1;	%persistencia choque cambial
// rhocstar,0.85;	%persistencia choque demanda externa.
// %rhom, beta_pdf, 0.85, 0.1; 		%persistencia choque monetario
// rhopid,0.85;		%mark-up choque no setor domestico
// rhopif,0.85;		%mark-up choque no setor importados
// rhopix,0.85;		%mark-up choque no setor exportados
// rhoc,0.85;		%mark-up choque no setor exportados

// // stderr epsm, 0.1;
// // stderr epsrstar, .1;
// // stderr epsveta, .1;
// // stderr epspid, .1;
// // stderr epspif, .1;
// // stderr epspix, .1;
// // stderr epscstar, .1;
// // stderr epsc, .1;
  
identification;
options_.gmhmaxlik.number=20000;
options_.gmhmaxlik.nscale=20000;
options_.gmhmaxlik.nclimb =0;
estimation(datafile=data,xls_sheet=data_77,xls_range=b1:j120,first_obs=1,mode_compute=6,mh_nblocks=3,mh_replic=1000000,lik_init=1,presample=0,mh_jscale=1);	

%impulse response

% shocks;
 //var epsveta; stderr 2.1/3;
 //var epscstar; stderr 0.2/3;
 //var epsm; stderr 0.45/3;
  //var epsrstar; stderr 0.15/3;
  //var epspid; stderr 0.9;
  //var epspif; stderr 0.9/3;
  //var epspix; stderr 1/2*0.9/3;
  //var epsc; stderr 1/2*0.6/3;
%end;

%stoch_simul(irf=10,nograph);