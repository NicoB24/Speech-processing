function var= var_cal(entradas,vector,media)
%Debo pasarle los vectores como fila!!
var=0;
for i=1:entradas
    var=var + (transpose((vector(i,:)-media)))*(vector(i,:)-media);
end
var=var/entradas;
end

