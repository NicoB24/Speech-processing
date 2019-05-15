%EM, aprendizaje no supervisado
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
%Las muestras totales de entrenamiento son
nro_entrenamiento=105;

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
%Creo vectores vacios para el LLH
LLH=[];
stop=0;
pos=0;

%Calculo las varianzas iniciales con las 5 muestras que uso para el bootstrap para
%cada letra
muestras_ai=Muestras_A([(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento],[1:Nro_formantes]);
muestras_oi=Muestras_O([(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento],[1:Nro_formantes]);
muestras_ui=Muestras_U([(Nro_muestras_entrenamiento-4):Nro_muestras_entrenamiento],[1:Nro_formantes]);
vari_a=var_calc(5, muestras_ai ,mediai_A);
vari_o=var_calc(5, muestras_oi ,mediai_O);
vari_u=var_calc(5, muestras_ui ,mediai_U);

%Calculo las varianzas iniciales con las muestras de entrenamientovari_a=var_calc(long_muestras,muestras,media_general)/2
% media_general=mean_calc(long_muestras,muestras);
% vari_a=var_calc(long_muestras,muestras,media_general)/3;
% vari_o=var_calc(long_muestras,muestras,media_general)/3;
% vari_u=var_calc(long_muestras,muestras,media_general)/3;

%Los pi_k son equiprobables al inicio
pi_a=1/3;
pi_o=1/3;
pi_u=1/3;

%Comienza la iteracion de EM
while (stop==0)
    pos=pos+1;
    gamma_veca=[];
    gamma_veco=[];
    gamma_vecu=[];
    for i=1:long_muestras
        gamma_a=mvnpdf(muestras(i,:),mediai_A,vari_a)*pi_a;
        gamma_o=mvnpdf(muestras(i,:),mediai_O,vari_o)*pi_o;
        gamma_u=mvnpdf(muestras(i,:),mediai_U,vari_u)*pi_u;
        denom = gamma_a+gamma_o+gamma_u;
        gamma_a = gamma_a/denom;
        gamma_o = gamma_o/denom;
        gamma_u = gamma_u/denom;
        gamma_veca=[gamma_veca; gamma_a];
        gamma_veco=[gamma_veco; gamma_o];
        gamma_vecu=[gamma_vecu; gamma_u]; 
    end
    
    %Recalculo las medias para cada clase
    mediai_A=mean_calcEM(long_muestras,muestras,gamma_veca);
    mediai_O=mean_calcEM(long_muestras,muestras,gamma_veco);
    mediai_U=mean_calcEM(long_muestras,muestras,gamma_vecu);
    %Recalculo la varianza para cada clase
    vari_a=var_calEM(long_muestras,muestras,mediai_A,gamma_veca);
    vari_o=var_calEM(long_muestras,muestras,mediai_O,gamma_veco);
    vari_u=var_calEM(long_muestras,muestras,mediai_U,gamma_vecu);
    %Recalculo los pi_k
    pi_a=pikEM(long_muestras,gamma_veca);
    pi_o=pikEM(long_muestras,gamma_veco);
    pi_u=pikEM(long_muestras,gamma_vecu);
    %Calculo el LLH
    px=0;
    for i=1:long_muestras
        pp=(mvnpdf(muestras(i,:),mediai_A,vari_a)*pi_a) + (mvnpdf(muestras(i,:),mediai_O,vari_o)*pi_o) + (mvnpdf(muestras(i,:),mediai_U,vari_u)*pi_u);
        lpp=log(pp);
        px=px+lpp;
    end
    LLH_new=px;    
    LLH=[LLH LLH_new];
    if length(LLH) > 1 && abs( LLH(end) - LLH(end-1) ) < 0.00001
        stop=1;
    end    
    %Grafico
    if(pos<=10)
        figure(pos+1)
        for i=1:long_muestras-1
            hold on
            plot(muestras(i,1),muestras(i,2),'o','color',[gamma_veca(i) gamma_vecu(i) gamma_veco(i)])
        end
        plot(mediai_A(:,1),mediai_A(:,2),'xr','linewidth',5,'markersize',10)
        plot(mediai_O(:,1),mediai_O(:,2),'xb','linewidth',5,'markersize',10)
        plot(mediai_U(:,1),mediai_U(:,2),'xg','linewidth',5,'markersize',10)
        title(sprintf('Muestras y medias de cada iteración con LLH de %d',LLH_new))
        hold off
    end
end


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
GKA=(-1/2)*log(det(vari_a))-(1/2)* transpose((Muestra_analisis-mediai_A)) * ((vari_a)^-1)*(Muestra_analisis-mediai_A) + log(pi_a);
GKO=(-1/2)*log(det(vari_o))-(1/2)* transpose((Muestra_analisis-mediai_O)) * ((vari_o)^-1)*(Muestra_analisis-mediai_O) + log(pi_o);
GKU=(-1/2)*log(det(vari_u))-(1/2)* transpose((Muestra_analisis-mediai_U)) * ((vari_u)^-1)*(Muestra_analisis-mediai_U) + log(pi_u);

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
figure(12)
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
elipseA = abs(chol(vari_a))*y';
elipse_corrida_A = elipseA' + transpose(mediai_A);
plot(elipse_corrida_A(:,1),elipse_corrida_A(:,2),'r');
%Elipse O
elipseO = abs(chol(vari_o))*y';
elipse_corrida_O = elipseO' + transpose(mediai_O);
plot(elipse_corrida_O(:,1),elipse_corrida_O(:,2),'b');
%Elipse U
elipseU = abs(chol(vari_u))*y';
elipse_corrida_U = elipseU' + transpose(mediai_U);
plot(elipse_corrida_U(:,1),elipse_corrida_U(:,2),'g');
legend('Media de A','Tests de A','Clasificaciones de A','Media de O','Tests de O','Clasificaciones de O','Media de U','Tests de U','Clasificaciones de U','Elipse de A','Elipse de O','Elipse de U','Location','southeast')
title('Muestras de test, muestras clasificadas y medias obtenidas de entrenamiento')
hold off