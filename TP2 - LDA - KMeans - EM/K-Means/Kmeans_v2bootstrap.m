%K-means, aprendizaje no supervisado
clear all 
close all
%Lectura de los archivos
Muestras_A=load('a.txt');
Muestras_O=load('o.txt');
Muestras_U=load('u.txt');
%Tomo parte de las muestras para entrenamiento.
Nro_muestras_entrenamiento=40;
%Numero de formantes a utilizar
Nro_formantes=2;
%Determino el numero de clases
clases=3;

%Inicializacion con bootstrap, tomo 5 muestras de las 40 que tengo para
%entrenamiento, primero las mezclo para tomar 5 aleatorias
Muestras_A=Muestras_A(randperm(length(Muestras_A)),:);
Muestras_O=Muestras_O(randperm(length(Muestras_O)),:);
Muestras_U=Muestras_U(randperm(length(Muestras_U)),:);
mediai_A=0;
mediai_O=0;
mediai_U=0;
for i=(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento
    mediai_A=mediai_A+Muestras_A(i,[1:Nro_formantes]);
end
mediai_A=mediai_A/5;
for i=(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento
    mediai_O=mediai_O+Muestras_O(i,[1:Nro_formantes]);
end
mediai_O=mediai_O/5;
for i=(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento
    mediai_U=mediai_U+ Muestras_U(i,[1:Nro_formantes]);
end
mediai_U=mediai_U/5;


%Armo un vector con todas las muestras no usadas para el bootstrap y lo
%mezclo
muestras=[Muestras_A([1:Nro_muestras_entrenamiento-5],[1:Nro_formantes]);Muestras_O([1:Nro_muestras_entrenamiento-5],[1:Nro_formantes]);Muestras_U([1:Nro_muestras_entrenamiento-5],[1:Nro_formantes])];
long_muestras=length(muestras);
perm=randperm(long_muestras);
muestras=muestras(perm,:);
%Vector donde voy a guardar las distorciones
D_vector=0;


%Comienzo a calcular las medias por iteración hasta que la distorción se
%estabilice
for p=1:6
    %Creo un vector de ceros para cada clase para ir concatenando
    claseA=zeros(1,Nro_formantes);
    claseO=zeros(1,Nro_formantes);
    claseU=zeros(1,Nro_formantes);
    %Primero calculo la distancia euclidea entre la media y la muestra, y luego
    %coloco a la muestra en la clase cuya distancia fue menor
    for i=1:long_muestras
        vect=muestras(i,:);
        z_A=(norm( vect - mediai_A ))^2;
        z_O=(norm( vect - mediai_O ))^2;
        z_U=(norm( vect - mediai_U ))^2;
        if(z_A < z_O) && (z_A < z_U)
            claseA=[claseA;vect];
        end
        if(z_O < z_A) && (z_O < z_U)
            claseO=[claseO;vect];
        end
        if(z_U < z_A) && (z_U < z_O)
            claseU=[claseU;vect];
        end
    end

    %Elimino el 0 del principio que agrego para iterar
    claseA=claseA([2:end],:);
    claseO=claseO([2:end],:);
    claseU=claseU([2:end],:); 
    %Calculo la cantidad de elementos de cada clase
    cant_A=length(claseA);
    cant_O=length(claseO);
    cant_U=length(claseU);

    %Actualizo los centroides
    mediai_A=mean_calc(cant_A,claseA);
    mediai_O=mean_calc(cant_O,claseO);
    mediai_U=mean_calc(cant_U,claseU);
    
    %Calculo distorcion para cada clase y la distorción total
    D_A=distor_calc(cant_A,claseA,mediai_A,clases);
    D_O=distor_calc(cant_O,claseO,mediai_O,clases);
    D_U=distor_calc(cant_U,claseU,mediai_U,clases);
    D=D_A+D_O+D_U;
    D_vector=[D_vector;D];
    
    figure(p)
    hold on
    plot(mediai_A(:,1),mediai_A(:,2),'xr','linewidth',5,'markersize',10)
    plot(claseA(:,1),claseA(:,2),'or')
    plot(mediai_O(:,1),mediai_O(:,2),'xb','linewidth',5,'markersize',10)
    plot(claseO(:,1),claseO(:,2),'ob')
    plot(mediai_U(:,1),mediai_U(:,2),'xg','linewidth',5,'markersize',10)
    plot(claseU(:,1),claseU(:,2),'og')
    legend('Media de a','Formantes a.txt','Media de o','Formantes o.txt','Media de u','Formantes u.txt','Location','southeast')
    title(sprintf('Formantes y medias en cada iteración con una distorcion de %d',D))
    hold off
end
%Tomo a partir del segundo componente ya que lo inicie en 0 para poder ir
%concatenando valores
D_vector=D_vector(2:end);
figure(p+1)
plot(D_vector)
title('Distorción')
ylabel('Valor')
xlabel('Iteracion')


%Calculo los pi_k
pi_a=cant_A/(3*(Nro_muestras_entrenamiento-5));
pi_o=cant_O/(3*(Nro_muestras_entrenamiento-5));
pi_u=cant_U/(3*(Nro_muestras_entrenamiento-5));
%Calculo la varianza
varianza_A=var_calc(cant_A, claseA,mediai_A);
varianza_O=var_calc(cant_O, claseO,mediai_O);
varianza_U=var_calc(cant_U, claseU,mediai_U);
%Armo un vector de muestras para probar el rendimiento.
Muestras_prueba=[Muestras_A([41:50],[1:Nro_formantes]);Muestras_O([41:50],[1:Nro_formantes]);Muestras_U([41:50],[1:Nro_formantes])];
Muestras_prueba_t=transpose(Muestras_prueba);


%Transpongo las medias para tenerlas como vectores columna
mediai_A=transpose(mediai_A);
mediai_O=transpose(mediai_O);
mediai_U=transpose(mediai_U);
%Comienzo a clasificar
clase_a_clasif=[];
clase_o_clasif=[];
clase_u_clasif=[];
i=1;
while(i<=30)
%Calculo las funciones g
Muestra_analisis=Muestras_prueba_t(:,i);
GKA=(-1/2)*log(det(varianza_A))-(1/2)* transpose((Muestra_analisis-mediai_A)) * ((varianza_A)^-1)*(Muestra_analisis-mediai_A) + log(pi_a);
GKO=(-1/2)*log(det(varianza_O))-(1/2)* transpose((Muestra_analisis-mediai_O)) * ((varianza_O)^-1)*(Muestra_analisis-mediai_O) + log(pi_o);
GKU=(-1/2)*log(det(varianza_U))-(1/2)* transpose((Muestra_analisis-mediai_U)) * ((varianza_U)^-1)*(Muestra_analisis-mediai_U) + log(pi_u);

%1 clase A - 2 clase O - 3 clase U

if(GKA>GKO)&&(GKA>GKU)
    respuesta=1;
    clase_a_clasif=[clase_a_clasif; transpose(Muestra_analisis)];
else if(GKO>GKA)&&(GKO>GKU)
    respuesta=2;
    clase_o_clasif=[clase_o_clasif; transpose(Muestra_analisis)];
    else
        respuesta=3;
        clase_u_clasif=[clase_u_clasif; transpose(Muestra_analisis)];
    end
end

Resultados(i)=respuesta;
i=i+1;
end

%Calculo los errores
errores=0;
for i=1:10
    if (Resultados(i)==2) ||(Resultados(i)==3)
        errores=errores+1;
    end
end
for i=11:20
    if (Resultados(i)==1) ||(Resultados(i)==3)
        errores=errores+1;
    end
end
for i=21:30
    if (Resultados(i)==1) ||(Resultados(i)==2)
        errores=errores+1;
    end
end
errores
porcentaje_error=errores/30



%Plot de las muestras de test, las muestras de test clasificas y las medias
figure(p+2)
hold on
plot(mediai_A(1,:),mediai_A(2,:),'xr','linewidth',5,'markersize',10)
plot(Muestras_A([41:50],1),Muestras_A([41:50],2),'.r')
plot(clase_a_clasif(:,1),clase_a_clasif(:,2),'or')
plot(mediai_O(1,:),mediai_O(2,:),'xb','linewidth',5,'markersize',10)
plot(Muestras_O([41:50],1),Muestras_O([41:50],2),'.b')
plot(clase_o_clasif(:,1),clase_o_clasif(:,2),'ob')
plot(mediai_U(1,:),mediai_U(2,:),'xg','linewidth',5,'markersize',10)
plot(Muestras_U([41:50],1),Muestras_U([41:50],2),'.g')
plot(clase_u_clasif(:,1),clase_u_clasif(:,2),'og')
%Elipses
t = 0:0.01:2*pi;
y =  [cos(t') sin(t')];
%Elipse A
elipseA = abs(chol(varianza_A))*y';
elipse_corrida_A = elipseA' + transpose(mediai_A);
plot(elipse_corrida_A(:,1),elipse_corrida_A(:,2),'r');
%Elipse O
elipseO = abs(chol(varianza_O))*y';
elipse_corrida_O = elipseO' + transpose(mediai_O);
plot(elipse_corrida_O(:,1),elipse_corrida_O(:,2),'b');
%Elipse U
elipseU = abs(chol(varianza_U))*y';
elipse_corrida_U = elipseU' + transpose(mediai_U);
plot(elipse_corrida_U(:,1),elipse_corrida_U(:,2),'g');
legend('Media de A','Tests de A','Clasificaciones de A','Media de O','Tests de O','Clasificaciones de O','Media de U','Tests de U','Clasificaciones de U','Elipse de A','Elipse de O','Elipse de U','Location','southeast')
title('Muestras de test, muestras clasificadas y medias obtenidas de entrenamiento')
hold off