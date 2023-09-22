function [IR]=ir_gensys(G1,impact,T,magn)
%G1 eh a matriz G1 do Gensys
%impact eh a matriz 'impact' do Gensys
%T eh o numero de periodos desejados
% magn eh o tamanho dos choques (opcional, padrao=1)

neq=size(G1,1); nchoques=size(impact,2);
IR=zeros(T,neq,nchoques);
G1=G1'; 
impact=impact';

if nargin<4
    magn=ones(nchoques,1);
end
    
for t=1:T
    for j=1:nchoques
        IR(t,:,j)=(magn(j)*impact(j,:))*G1^(t-1);
    end
end