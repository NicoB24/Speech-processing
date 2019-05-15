function mean = mean_calc(entradas,vector)
mean=0;
for i=1:entradas
    mean=mean+vector(i,:);
end
mean=mean/entradas;
end
