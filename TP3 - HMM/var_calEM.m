function var= var_calEM(entradas,vector,media,gamma)
%Debo pasarle los vectores como fila!!
var=0;
sum_gamma=0;

for i=1:entradas
    var=var + gamma(i)*(transpose((vector(i,:)-media)))*(vector(i,:)-media);
end

for i=1:entradas
    sum_gamma=sum_gamma+gamma(i);
end

var=var/sum_gamma;

end

