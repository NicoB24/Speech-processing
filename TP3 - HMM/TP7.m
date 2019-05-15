%%
%Ejercicio 1

inic_hmm
%EM
means4 = hmm4.means;
vars4= hmm4.vars;
means = hmm4.means;
vars = hmm4.vars;
%Matriz de transicion inicial
trans = [0   1   0   0   0  ;...
         0   0.5 0.5 0   0  ;...
         0   0   0.5 0.5 0  ;...
         0   0   0   0.5 0.5;...
         0   0   0   0   1  ];
disp('Ejercicio 1')
while true
    [x,stateSeq] = genhmm(hmm4.means,hmm4.vars,hmm4.trans);
    
    %Me aseguro que x tenga 100 elementos
    if length(x) < 100
        continue
    end
    
    %Me aseguro que hayan 15 elementos en cada clase
    for z = 2:length(hmm4.means)-1
        if sum(stateSeq == z) < 15
            continue
        end
    end
    
    break
    
end

%Calculo la media y varianza inicial
longitud=length(x);
media_inicial=mean_calc(longitud,x);
varianza_inicial=var_calc(longitud,x,media_inicial);
%Armo los arrays que luego utilizo para calcular los parametros en logfwd
media_inicial=transpose(media_inicial);
means(2:end-1) = {media_inicial};    
vars(2:end-1) = {varianza_inicial};    


likelihood = [];
stop=0;
iteracion=1;
while(stop==0)

    %Calculo los parametros de Markov con logfwd
    [logProb_alpha,alpha_matrix,logProb_beta,beta_matrix,gamma_matrix,xi_matrix] = logfwd2(x,means,vars,trans);
    gamma_matrix = exp(gamma_matrix)'; %Hago la exponencial a la gamma y tranpongo para que me queden como columnas
    
    
    %Recalculo la media y la varianza
    for k = 2:length(means)-1        
        longitud=length(x);
        gamma=gamma_matrix(:,k);
        
        %Media
        media=mean_calcEM(longitud,x,gamma);
        means{k}=transpose(media);
        
        %Varianza
        varianza=var_calEM(longitud,x,media,gamma);
        if (rcond(varianza)>0.000000000001)
            vars{k}=varianza;
        end
    end

    
    %Calculo la matriz de transicion
    xi_normal = exp(xi_matrix);
    for j = 1:size(xi_normal,1)
        for k = 1:size(xi_normal,2)           
            numerador = sum(xi_normal(j,k,2:end));
            denominador = sum(sum(xi_normal(j,:,2:end)));
            trans(j+1,k+1) = numerador/denominador;          
        end
    end
    
    trans(end-1,end) = abs(1 - sum(trans(end-1,1:end-1)));
    
    %Dibujo las primeras 25 iteraciones (a veces itera demasiadas veces y
    %se cuelga)
    if(iteracion<15)
        figure
        hold on
        colores = 'rgb';
        %Puntos
        for t = 1:size(x,1)
            aux = gamma_matrix(t,2:end)*1000;
            aux = floor(aux)/1000;
            plot(x(t,1),x(t,2), 'o', 'color', aux)
        end
        %Elipses
        for k = 2:length(means)-1
            t = 0:0.01:2*pi;
            y =  [cos(t') sin(t')];
            elipse = abs(chol(vars{k}))*y';
            elipse_corrida = elipse' + transpose(means{k});
            plot(elipse_corrida(:,1),elipse_corrida(:,2),colores(k-1));
            plot(means{k}(1), means{k}(2), ['+' colores(k-1)]);
        end
    end
    title('Muestras clasificadas, medias obtenidas y curvas de nivel')
    
    %Analizo si la probabilidad se mantiene, sino, itero de vuelta
    likelihood=[likelihood logProb_alpha];
    if length(likelihood) > 1 && abs( likelihood(end) - likelihood(end-1) ) < 0.0000001
        stop=1;
    end    
    iteracion=iteracion+1;
end
hold off

figure
%Grafico del modelo original
plotseq2(x,stateSeq)
title('Modelo original')
%Matriz de transicion original y estimada
disp('La matriz de transicion original y la estimada son:')
trans_original=hmm4.trans
trans

%Grafico de las medias y varianzas originales y estimadas
figure
hold on
colores = 'rgb';
 for k = 2:length(means)-1
    t = 0:0.01:2*pi;
    y =  [cos(t') sin(t')];
   
    elipse_original = abs(chol(vars4{k}))*y';
    elipse_corrida_original = elipse_original' + transpose(means4{k});
    plot(elipse_corrida_original(:,1),elipse_corrida_original(:,2),colores(k-1));
    
    elipse = abs(chol(vars{k}))*y';
    elipse_corrida = elipse' + transpose(means{k});
    plot(elipse_corrida(:,1),elipse_corrida(:,2),colores(k-1));
    
    plot(means4{k}(1), means4{k}(2), ['+' colores(k-1)]);
    plot(means{k}(1), means{k}(2), ['+' colores(k-1)]);
 end
title('Medias y varianzas originales y estimadas')

%%
%Ejercicio 2
disp('Ejercicio 2')
%Parametros de hmm4
mean4 = hmm4.means;
var4 = hmm4.vars;
trans4 = hmm4.trans;
%Parametros de hmm6
mean6 = hmm6.means;
var6 = hmm6.vars;
trans6 = hmm6.trans;

%Armo la matriz de transicion del modelo
trans(2:4,2:4)=trans4(2:4,2:4);
trans(5:7,5:7)=trans6(2:4,2:4);
%Igual probabilidad de empezar en uno u otro
trans(1,2)=0.5;
trans(1,5)=0.5;
%Arreglo los otros términos
trans(4,2)=0.05/3;
trans(4,4)=0.95;
trans(4,5)=0.05/3;
trans(4,8)=0.05/3;
trans(7,2)=0.05/3;
trans(7,5)=0.05/3;
trans(7,7)=0.95;
trans(7,8)=0.05/3;
trans(8,8)=1;

disp('La matriz de transicion es')
trans

%Armo la media y varianza del modelo
l=length(mean4);
mean=mean4;
var=var4;
%Armo la media
for i=l:(l+5-2)
    mean{i}=mean6{i-l+2};
end
%Armo la varianza
for i=l:(l+5-2)
    var{i}=var6{i-l+2};
end

%Armo una secuencia de observaciones y la calculo por Viterbi
[x , stateSeq] = genhmm(mean,var,trans);
[stateSeq_viterbi , logProb] = logvit(x,mean,var,trans);
if(isequal(stateSeq,stateSeq_viterbi))
    disp('La secuencia obtenida por Viterbi y la que realmente se produjo son iguales')
end

%Comparo la probabilidad calculada por Viterbi con la probabilidad real 
[logProb_alpha,alpha_matrix,logProb_beta,beta_matrix,gamma_matrix,xi_matrix] = logfwd2(x,mean,var,trans);
disp('La probabilidad obtenida utilizando logfwd y Viterbi son respecticamente:')
logProb_alpha
logProb

%Analizo la secuencia de modelos que se produjo
sequence=[];
for i=2:(length(stateSeq)-1)
    if( stateSeq(i)>=2 && stateSeq(i)<=4)
        sequence=[sequence 4];
    end
    if( stateSeq(i)>=5 && stateSeq(i)<=7)
        sequence=[sequence 6];
    end
end

disp('La secuencia de modelos ocurridas es:')
sequence
