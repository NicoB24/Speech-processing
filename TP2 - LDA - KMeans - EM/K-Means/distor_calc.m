function D = distor_calc( entradas,vector,media, clases )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
D=0;
for i=1:entradas
    D=D+((norm(vector(i,:)-media)))^2;
end
D=D/entradas;
end

