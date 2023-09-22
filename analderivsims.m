function mod=anal_deriv(mod);
%Calcula as derivadas analiticas do modelo LINEARIZADO para computar as
%matrizes do GENSYS
%Inputs: estrutura mod contendo o modelo com as equacoes escritas no
%formato g0*y(t+1)=g1y(t) + PSI eps + PI eta
%f eh o sistema de eqs.
%eps sao choques exogenos
%eta sao erros expectacionais

y = mod.Y;
yf = mod.Yf;
eps = mod.eps;
et = mod.et;
f = mod.f;
%Matrizes g0 e g1 PI e PSI do GENSYS

g1 = jacobian(f,y);
g0 = jacobian(f,yf);
PI = jacobian(f,et);
PSI = jacobian(f,eps);

mod.g1 = g1;
mod.g0 = g0;
mod.PI = PI;
mod.PSI = PSI;
end
