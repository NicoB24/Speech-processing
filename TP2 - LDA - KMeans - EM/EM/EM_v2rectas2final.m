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

%Inicializacion aleatoria
%Mezclo las muestras
Muestras_A=Muestras_A(randperm(length(Muestras_A)),:);
Muestras_O=Muestras_O(randperm(length(Muestras_O)),:);
Muestras_U=Muestras_U(randperm(length(Muestras_U)),:);

%Calculo una media general
muestras_totales=[Muestras_A;Muestras_O;Muestras_U];
%Me quedo con las dos primeras columnas de las muestras
muestras_totales=muestras_totales(:,[1:Nro_formantes]);
long_muestras_totales=length(muestras_totales);
media_total=mean_calc(long_muestras_totales,muestras_totales);
%Me quedo con las dos primeras componentes de la media
media_total=media_total(:,[1:Nro_formantes]);
%Me quedo con las 40 de entrenamiento
muestras_ai=Muestras_A([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
muestras_oi=Muestras_O([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
muestras_ui=Muestras_U([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
muestras_i=[muestras_ai;muestras_oi;muestras_ui];
long_i=length(muestras_i);

%Comienzo a clasificar
claseai=[];
claseoi=[];
claseui=[];
while ( isempty(claseai) ) || ( isempty(claseoi) ) || ( isempty(claseui) )
%Defino las clases
claseai=zeros(1,Nro_formantes);
claseoi=zeros(1,Nro_formantes);
claseui=zeros(1,Nro_formantes);
%Armo 3 espacios aleatorios separados 120 grados
angle=360*rand;
angulo1=angle+120;
angulo2=angulo1+120;
if(angulo1>360)
    angulo1=angulo1-360;
end
if(angulo2>360)
    angulo2=angulo2-360;
end
angulos_i=[angle;angulo1;angulo2];
long_angulos_i=length(angulos_i);
max_angle=max(angulos_i);
min_angle=min(angulos_i);
for m=1:long_angulos_i
    valor_angulo=angulos_i(m);
    if (valor_angulo~=max_angle) && (valor_angulo~=min_angle)
        mid_angle=valor_angulo;
    end
end
    
for j=1:long_i
      resta=muestras_i(j,:)-media_total;
      angulo_r=myatan(resta(1,2),resta(1,1));
      if (angulo_r>=min_angle) && (angulo_r<=mid_angle)
        claseai=[claseai;muestras_i(j,:)];
      end
      if (angulo_r>mid_angle) && (angulo_r<=max_angle)
        claseoi=[claseoi;muestras_i(j,:)];
      end
      if (angulo_r>max_angle) || (angulo_r<min_angle)
        claseui=[claseui;muestras_i(j,:)];
      end
end
claseai=claseai([2:end],:);
claseoi=claseoi([2:end],:);
claseui=claseui([2:end],:);
[longa,xa]=size(claseai);
[longo,xo]=size(claseoi);
[longu,xu]=size(claseui);
%El If siguiente lo uso para evitar que una clase quede con un elemento
%solo
if (longa==1) || (longo==1) || (longu==1)
    claseai=[];
    claseoi=[];
    claseui=[];
end
end
%Calculo las medias iniciales 
mediai_A=mean_calc(longa,claseai);
mediai_O=mean_calc(longo,claseoi);
mediai_U=mean_calc(longu,claseui);

%Graficos
puntos=linspace(1,2000);
mediax=media_total(:,1)*ones(1,100);
mediay=media_total(:,2)*ones(1,100);
figure(1)
hold on
plot(puntos*cosd(angle)+mediax,puntos*sind(angle) +mediay,'k')
plot(puntos*cosd(angulo1)+mediax,puntos*sind(angulo1) +mediay,'k')
plot(puntos*cosd(angulo2)+mediax,puntos*sind(angulo2) +mediay,'k')
plot(media_total(:,1),media_total(:,2),'xk','linewidth',5,'markersize',10)
plot(mediai_A(:,1),mediai_A(:,2),'xr','linewidth',5,'markersize',10)
plot(claseai(:,1),claseai(:,2),'or')
plot(mediai_O(:,1),mediai_O(:,2),'xb','linewidth',5,'markersize',10)
plot(claseoi(:,1),claseoi(:,2),'ob')
plot(mediai_U(:,1),mediai_U(:,2),'xg','linewidth',5,'markersize',10)
plot(claseui(:,1),claseui(:,2),'og')
legend('Recta angle','Recta angulo1','Recta angulo2','Media total','Media de a','Formantes a.txt','Media de o','Formantes o.txt','Media de u','Formantes u.txt','Location','northeast');
title('Formantes tomados al azar y medias iniciales')
hold off


%Armo un vector con todas las muestras de entrenamiento y las mezclo
muestras=[Muestras_A([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);Muestras_O([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);Muestras_U([1:Nro_muestras_entrenamiento],[1:Nro_formantes])];
long_muestras=length(muestras);
perm=randperm(long_muestras);
muestras=muestras(perm,:);
%Creo vectores vacios para el LLH
LLH=-100;
stop=0;
pos=1;

%Calculo las varianzas con las muestras de entrenamiento
vari_a=var_calc(longa, claseai ,mediai_A);
vari_o=var_calc(longo, claseoi ,mediai_O);
vari_u=var_calc(longu, claseui ,mediai_U);
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


%Calculo las medias de las muestras de testeo
media_A_test=mean_calc(10,Muestras_A([41:50],[1:Nro_formantes]));
media_O_test=mean_calc(10,Muestras_O([41:50],[1:Nro_formantes]));
media_U_test=mean_calc(10,Muestras_U([41:50],[1:Nro_formantes]));
%Tengo 3 medias finales
media_Af=0;
media_Of=0;
media_Uf=0;
%Me fijo a que media se parece mas cada uno
%A
resta1=(norm( media_A_test - mediai_A ))^2;
resta2=(norm( media_A_test - mediai_O ))^2;
resta3=(norm( media_A_test - mediai_U ))^2;
[minA,iminA] = min([resta1 resta2 resta3]);
if (iminA == 1)
    media_Af=mediai_A;
    varianza_Af=vari_a;
else if (iminA == 2)
    media_Af=mediai_O;
    varianza_Af=vari_o;
    else
        media_Af=mediai_U;
        varianza_Af=vari_u;
    end
end

%O
resta1=(norm( media_O_test - mediai_A ))^2;
resta2=(norm( media_O_test - mediai_O ))^2;
resta3=(norm( media_O_test - mediai_U ))^2;
[minA,iminA] = min([resta1 resta2 resta3]);
if (iminA == 1)
    media_Of=mediai_A;
    varianza_Of=vari_a;
else if (iminA == 2)
    media_Of=mediai_O;
    varianza_Of=vari_o;
    else
        media_Of=mediai_U;
        varianza_Of=vari_u;
    end
end

%U
resta1=(norm( media_U_test - mediai_A ))^2;
resta2=(norm( media_U_test - mediai_O ))^2;
resta3=(norm( media_U_test - mediai_U ))^2;
[minA,iminA] = min([resta1 resta2 resta3]);
if (iminA == 1)
    media_Uf=mediai_A;
    varianza_Uf=vari_a;
else if (iminA == 2)
    media_Uf=mediai_O;
    varianza_Uf=vari_o;
    else
        media_Uf=mediai_U;
        varianza_Uf=vari_u;
    end
end



%Transpongo las medias para tenerlas como vectores columna
media_Af=transpose(media_Af);
media_Of=transpose(media_Of);
media_Uf=transpose(media_Uf);

%Comienzo a clasificar
clase_a_clasif=[];
clase_o_clasif=[];
clase_u_clasif=[];
i=1;
while(i<=30)
%Calculo las funciones g
Muestra_analisis=Muestras_prueba_t(:,i);
GKA=(-1/2)*log(det(varianza_Af))-(1/2)* transpose((Muestra_analisis-media_Af)) * ((varianza_Af)^-1)*(Muestra_analisis-media_Af) + log(pi_a);
GKO=(-1/2)*log(det(varianza_Of))-(1/2)* transpose((Muestra_analisis-media_Of)) * ((varianza_Of)^-1)*(Muestra_analisis-media_Of) + log(pi_o);
GKU=(-1/2)*log(det(varianza_Uf))-(1/2)* transpose((Muestra_analisis-media_Uf)) * ((varianza_Uf)^-1)*(Muestra_analisis-media_Uf) + log(pi_u);

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
tasa_de_error=errores/30


%Plot de las muestras de test, las muestras de test clasificas y las medias
figure(12)
hold on
plot(media_Af(1,:),media_Af(2,:),'xr','linewidth',5,'markersize',10)
plot(Muestras_A([41:50],1),Muestras_A([41:50],2),'.r')
plot(clase_a_clasif(:,1),clase_a_clasif(:,2),'or')
plot(media_Of(1,:),media_Of(2,:),'xb','linewidth',5,'markersize',10)
plot(Muestras_O([41:50],1),Muestras_O([41:50],2),'.b')
plot(clase_o_clasif(:,1),clase_o_clasif(:,2),'ob')
plot(media_Uf(1,:),media_Uf(2,:),'xg','linewidth',5,'markersize',10)
plot(Muestras_U([41:50],1),Muestras_U([41:50],2),'.g')
plot(clase_u_clasif(:,1),clase_u_clasif(:,2),'og')
% %Elipses
t = 0:0.01:2*pi;
y =  [cos(t') sin(t')];
%Elipse A
elipseA = abs(chol(varianza_Af))*y';
elipse_corrida_A = elipseA' + transpose(media_Af);
plot(elipse_corrida_A(:,1),elipse_corrida_A(:,2),'r');
%Elipse O
elipseO = abs(chol(varianza_Of))*y';
elipse_corrida_O = elipseO' + transpose(media_Of);
plot(elipse_corrida_O(:,1),elipse_corrida_O(:,2),'b');
%Elipse U
elipseU = abs(chol(varianza_Uf))*y';
elipse_corrida_U = elipseU' + transpose(media_Uf);
plot(elipse_corrida_U(:,1),elipse_corrida_U(:,2),'g');
legend('Media de A','Tests de A','Clasificaciones de A','Media de O','Tests de O','Clasificaciones de O','Media de U','Tests de U','Clasificaciones de U','Elipse de A','Elipse de O','Elipse de U','Location','southeast')
title('Muestras de test, muestras clasificadas y medias obtenidas de entrenamiento')
hold off
