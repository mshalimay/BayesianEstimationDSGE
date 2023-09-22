function [param,set] = parameters()
set.bet=0.99;
set.alphay=0;
set.alphac=0.3;
set.v=.85;
set.rrbar=1/set.bet;
set.phibar=1-set.v*set.bet;
set.omeg=2;
set.rhor=0.8;
set.phiy=0.125;
set.rhoc=0.9;
set.rhopif=0.9;

%% Parametros estimados
%chutes iniciais sao as medias das priors (ver priorarray para saber as dist)

%parametros de rigidez
param.xid=0.675; %domestic rigidity
param.xif=0.5;    %import rigidity
param.xix=0.5;    %export rigidity

%Parametros de indexacao
param.gammad=0.5;	%domestic indexation
param.gammaf=0.5;		%import
param.gammax=0.5;		%export

%Parametros comportamentais
param.sigm=1.5;   %elast. intertemp. subst.
param.eta=1.5;		%elast. subst. import x domestic goods

%Regra de taylor
param.phi=1.7;		%inflation response TR
param.rdpi=0.3;  	%inflation growth response TR
param.phirer=0; 	%real exchange rate response TR
param.rdy=0.0625;  %gdp growth response TR

%Autoregressivos dos choques
param.rhoveta=0.9;	  %productivity shock persistence
param.rhorstar=0.9; %exchange rate shock persistence
param.rhocstar=0.9;	  %external demand
param.rhom=0.9; 	%monetary shock
param.rhopix=0.9;		%mark-up export sector

%Desvios padroes dos choques

%inv gamma mean x (mean = 3*mode=scale) and 2 degree freedom (shape) (or sd = infinity)
param.sdm=0.45/3; 		 %moetary shock 
param.sdrstar=0.15/3;   %exchange rate shock 
param.sdvet=2.1/3;     %product. shock
param.sdpif=0.9/3;     %mark-up import
param.sdpix=0.9/3;     %mark-up export
param.sdcstar=0.6/3;   %foreign demand shock
param.sdc=0.6/3;       %preference shock

end
