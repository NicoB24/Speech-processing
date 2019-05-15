function pik= pikEM(entradas,gamma)
%Debo pasarle los vectores como fila!!
pik=0;

for i=1:entradas
    pik=pik+gamma(i);
end

pik=pik/entradas;

end