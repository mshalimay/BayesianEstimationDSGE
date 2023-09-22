function model_func(model)
f = fopen('model_prog.m', 'w');

%chama os parametros do modelo
param_list = char(model.PARAM);
param_list = param_list(10:end-3);
str = ['function [g0 g1 PSI PI] = model_prog(param, set)'];
fprintf(f, '%s\n', str);
fprintf(f, '%s\n', '');

%parametro para string
str_sv{1} = ['function ' str(1:end-1)];
str = '%parametro para string';
str_sv{2}= str;
fprintf(f, '%s\n', str);
for j = 1:length(model.PARAM)
    str = [char(model.PARAM(j)) ' = param(' num2str(j), ');'];
    str_sv{j+2} = str;
    fprintf(f, '%s\n',str);
end
fprintf(f, '%s\n', '');


%set para string
str = '%set para string';
str_sv{j+3} = str;
fprintf(f, '%s\n', str);
for j = 1:length(model.SET)
    str = [char(model.SET(j)) ' = set(' num2str(j), ');'];
    str_sv{j+length(model.PARAM)+4} = str;
    fprintf(f, '%s\n',str);
end
fprintf(f, '%s\n', '');

%***************************
%g0
%***************************
str = ['g0 = ' symmat_print(model.g0) ';'];
fprintf(f, '%s\n', str);

%***************************
%g1
%***************************
str = ['g1 = ' symmat_print(model.g1) ';'];
fprintf(f, '%s\n', str);

%***************************
%PSI
%***************************
str = ['PSI = ' symmat_print(model.PSI) ';'];
fprintf(f, '%s\n', str);

%***************************
%PI
%***************************
str = ['PI = ' symmat_print(model.PI) ';'];
fprintf(f, '%s\n', str);


fclose(f);

%*******************************************************
% Transforma em matriz do matlab
%*******************************************************
function str = symmat_print(x)

str = char(x);
str = str(8:end-1);

row_idx = findstr(str, '],');
for j = 1:length(row_idx)
    str(row_idx(j):row_idx(j)+1)='];';
end

%*******************************************************
%troca old por new_str
%*******************************************************
function update_param(old, new_str)

new = [old, '_tmp'];

f_old = fopen([old, '.m'], 'r');
f_new = fopen(new, 'w');

l = '';
j = 1;
while ~strcmp('%START_ADD', l) && j < 300
  l = fgetl(f_old);
  fprintf(f_new,'%s\n',l);
  j = j+1;
end
if j == 300
    error('Too many lines have passed.');
end

for j = 1:length(new_str)
  fprintf(f_new, '%s\n', new_str{j});
end

j = 1;
while ~strcmp('%END_ADD', l) && j < 300
    l = fgetl(f_old);
    j = j+1;
end
if j == 300
    error('Too many lines have passed.');
end
fprintf(f_new,'%s',l);

j = 1;
while ischar(l) && j < 500
  fprintf(f_new, '\n',''); 
  l = fgetl(f_old);
  fprintf(f_new,'%s',l); 
  j = j+1;
end
if j == 500
    error('Too many lines have passed.');
end

fclose(f_old);
fclose(f_new);
eval(['!mv ', new ' ' old '.m']);

