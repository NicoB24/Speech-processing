%LDA, entrenamiento supervisado.

%Lectura de los archivos
Muestras_A=load('a.txt');
Muestras_O=load('o.txt');
Muestras_U=load('u.txt');
%Tomo parte de las muestras para entrenamiento.
Nro_muestras_entrenamiento=40;
%Numero de formantes a utilizar
Nro_formantes=2;

%Selecciono en base a las muestras y formantes, las muestras de
%entrenamiento
Muestras_entrenamiento_A=Muestras_A([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
Muestras_entrenamiento_O=Muestras_O([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
Muestras_entrenamiento_U=Muestras_U([1:Nro_muestras_entrenamiento],[1:Nro_formantes]);
%Hago el grafico de los primeros 2 formantes para los tres tipos.
figure(1)
hold on
plot(Muestras_entrenamiento_A(:,1),Muestras_entrenamiento_A(:,2),'ro');
plot(Muestras_entrenamiento_O(:,1),Muestras_entrenamiento_O(:,2),'go');
plot(Muestras_entrenamiento_U(:,1),Muestras_entrenamiento_U(:,2),'bo');



%Calculo las medias de las muestras de entrenamiento para cada tipo.
%media_A
media_A=mean_calc(Nro_muestras_entrenamiento, Muestras_A(:,[1:Nro_formantes]));
%media_O
media_O=mean_calc(Nro_muestras_entrenamiento, Muestras_O(:,[1:Nro_formantes]));
%media_U
media_U=mean_calc(Nro_muestras_entrenamiento, Muestras_U(:,[1:Nro_formantes]));


%Calculo la varianza
%Utilizo el transpuesto al reves ya que los vectores en las formulas se piensan
%como vector columna
%Varianza_A
varianza_A=var_calc(Nro_muestras_entrenamiento,Muestras_A(:,[1:Nro_formantes]),media_A);
%Varianza_O
varianza_O=var_calc(Nro_muestras_entrenamiento,Muestras_O(:,[1:Nro_formantes]),media_O);
%Varianza_U
varianza_U=var_calc(Nro_muestras_entrenamiento,Muestras_U(:,[1:Nro_formantes]),media_U);
%La varianza a utilizar en LDA es un promedio ponderado de las 3
varianza_ML=(Nro_muestras_entrenamiento*varianza_A + Nro_muestras_entrenamiento*varianza_O + Nro_muestras_entrenamiento*varianza_U)/(3*Nro_muestras_entrenamiento);



%Para verifica a que categoria o clase pertenece uso la funcion g
%Elijo la muestras para analizar
%Muestra_analisis=transpose(Muestras_U(45,[1:Nro_formantes]));

%Armo un vector de muestras para probar el rendimiento.
Muestras_prueba=[Muestras_A([41:50],[1:Nro_formantes]);Muestras_O([41:50],[1:Nro_formantes]);Muestras_U([41:50],[1:Nro_formantes])];
Muestras_prueba_t=transpose(Muestras_prueba);

i=1;
while(i<=30)

Muestra_analisis=Muestras_prueba_t(:,i);    
XOk_A=media_A*((varianza_ML)^-1)*transpose(media_A);
WOk_A=transpose(((varianza_ML)^-1)*transpose(media_A));
GKA=(-1/2)*XOk_A+WOk_A*Muestra_analisis+log(1/3);

XOk_O=media_O*((varianza_ML)^-1)*transpose(media_O);
WOk_O=transpose(((varianza_ML)^-1)*transpose(media_O));
GKO=(-1/2)*XOk_O+WOk_O*Muestra_analisis+log(1/3);

XOk_U=media_U*((varianza_ML)^-1)*transpose(media_U);
WOk_U=transpose(((varianza_ML)^-1)*transpose(media_U));
GKU=(-1/2)*XOk_U+WOk_U*Muestra_analisis+log(1/3);

%1 clase A 2 clase O 3 clase U

if(GKA>GKO)&&(GKA>GKU)
    respuesta=1;
end
if(GKO>GKA)&&(GKO>GKA)
    respuesta=2;
end
if(GKU>GKO) && (GKU>GKA)
    respuesta=3;
end
Resultados(i)=respuesta;
i=i+1;
end


%Elipses
t = 0:0.01:2*pi;
y =  [cos(t') sin(t')];
%Elipse A
elipseA = abs(chol(varianza_A))*y';
elipse_corrida_A = elipseA' + media_A;
plot(elipse_corrida_A(:,1),elipse_corrida_A(:,2),'r');
%Elipse O
elipseO = abs(chol(varianza_O))*y';
elipse_corrida_O = elipseO' + media_O;
plot(elipse_corrida_O(:,1),elipse_corrida_O(:,2),'g');
%Elipse U
elipseU = abs(chol(varianza_U))*y';
elipse_corrida_U = elipseU' + media_U;
plot(elipse_corrida_U(:,1),elipse_corrida_U(:,2),'b');
legend('Formantes a.txt','Formantes o.txt','Formantes u.txt','Elipse de A','Elipse de O','Elipse de U');

