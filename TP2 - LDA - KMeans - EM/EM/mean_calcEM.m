function mean = mean_calcEM(entradas,vector,gamma)
mean=0;
sum_gamma=0;

for i=1:entradas
    mean=mean+vector(i,:)*gamma(i);
end

for i=1:entradas
    sum_gamma=sum_gamma+gamma(i);
end

mean=mean/sum_gamma;

end
