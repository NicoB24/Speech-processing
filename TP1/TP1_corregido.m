%%
%ITEMS DE VERIFICACION (items 1 y 2)

[audio,Fs]=audioread('fantasia.wav',[14000,14000+0.025*16000]); %Señal s
autocorrelacion=xcorr(audio);
%plot(autocorrelacion)
[maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo 19 muestras mas
matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
b=[0 ;(matriz_ro^-1)*vector_ro]; %Calculo los b=(matriz_p^-1)*p
estimacion=filter(b,1,audio); %Calculo mi señal s_sombero
error=audio-estimacion; %s - s_sombrero

figure(1)
hold on
plot(audio,'r')
plot(estimacion,'b')
title('Fragmento de audio original y su estimacion')
xlabel('Número de muestra')
ylabel('Amplitud')
legend('Señal de audio','Señal estimada','Location','southwest')
hold off

figure(2)
plot(error)
title('Señal del error')
xlabel('Número de muestra')
ylabel('Amplitud')


b_g=(matriz_ro^-1)*vector_ro; %b sin el corrimiento que uso para graficar
g=sqrt(autocorrelacion(max_pos)-(transpose(b_g))*vector_ro);
[h,w]=freqz(g,[1;((-1)*b_g)],1024,'whole');
L=1024; %Cantidad de puntos
f = Fs*(0:L-1)/L; %Eje de freucencias
fft_audio=fft(audio,1024);
fft_estimacion=fft(estimacion,1024);
figure(3)
hold on
plot(f,abs(h))
plot(f,abs(fft_audio),'r')
title('FFT de la señal de audio original y la envolvente')
xlabel('Frecuencia')
ylabel('Amplitud')
legend('Envolvente','FFT de la señal de audio','Location','northeast')
hold off


%%
%SUPERPOSICION DE ENVOLVENTES (item 3)

%Primera vocal a
[audio,Fs]=audioread('fantasia.wav',[11100,11500]); %Señal s
autocorrelacion=xcorr(audio);
%plot(autocorrelacion)
[maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo 19 muestras mas
matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
b=[0 ;(matriz_ro^-1)*vector_ro]; %Calculo los b=(matriz_p^-1)*p
estimacion=filter(b,1,audio); %Calculo mi aseñal s_sombero
error=audio-estimacion; %s - s_sombrero
b_g=(matriz_ro^-1)*vector_ro; %b sin el corrimiento que uso para graficar
g=sqrt(autocorrelacion(max_pos)-(transpose(b_g))*vector_ro);
[h1,w]=freqz(g,[1;((-1)*b_g)],1024,'whole');

%Segunda vocal a
[audio,Fs]=audioread('fantasia.wav',[14000,14400]); %Señal s
autocorrelacion=xcorr(audio);
%plot(autocorrelacion)
[maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo 19 muestras mas
matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
b=[0 ;(matriz_ro^-1)*vector_ro]; %Calculo los b=(matriz_p^-1)*p
estimacion=filter(b,1,audio); %Calculo mi aseñal s_sombero
error=audio-estimacion; %s - s_sombrero
b_g=(matriz_ro^-1)*vector_ro; %b sin el corrimiento que uso para graficar
g=sqrt(autocorrelacion(max_pos)-(transpose(b_g))*vector_ro);
[h2,w]=freqz(g,[1;((-1)*b_g)],1024,'whole');

%Tercera vocal i
[audio,Fs]=audioread('fantasia.wav',[18700,20100]); %Señal s
autocorrelacion=xcorr(audio);
%plot(autocorrelacion)
[maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo 19 muestras mas
matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
b=[0 ;(matriz_ro^-1)*vector_ro]; %Calculo los b=(matriz_p^-1)*p
estimacion=filter(b,1,audio); %Calculo mi aseñal s_sombero
error=audio-estimacion; %s - s_sombrero
b_g=(matriz_ro^-1)*vector_ro; %b sin el corrimiento que uso para graficar
g=sqrt(autocorrelacion(max_pos)-(transpose(b_g))*vector_ro);
[h3,w]=freqz(g,[1;((-1)*b_g)],1024,'whole');

%Cuarta vocal a
[audio,Fs]=audioread('fantasia.wav',[21400,21800]); %Señal s
autocorrelacion=xcorr(audio);
%plot(autocorrelacion)
[maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo 19 muestras mas
matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
b=[0 ;(matriz_ro^-1)*vector_ro]; %Calculo los b=(matriz_p^-1)*p
estimacion=filter(b,1,audio); %Calculo mi aseñal s_sombero
error=audio-estimacion; %s - s_sombrero
b_g=(matriz_ro^-1)*vector_ro; %b sin el corrimiento que uso para graficar
g=sqrt(autocorrelacion(max_pos)-(transpose(b_g))*vector_ro);
[h4,w]=freqz(g,[1;((-1)*b_g)],1024,'whole');

L=1024; %Cantidad de puntos
f = Fs*(0:L-1)/L; %Eje de freucencias
figure(4)
hold on
plot(f,abs(h1))
plot(f,abs(h2),'r')
plot(f,abs(h3),'g')
plot(f,abs(h4),'y')
title('Envolventes correspondistes a cada señal')
xlabel('Frecuencia')
ylabel('Amplitud')
legend('Primera vocal a','Segunda vocal a', 'Tercera vocal i', 'Cuarta vocal a', 'Location','northeast')
hold off

%%
%SUPERFICIE COMPLETA DE LAS ENVOLVENTES Y ESPECTROGRAMA (item 4)

i=1; %Posición en la columna de H
index=1; %Posición de lectura del audio
[y,Fs] = audioread('fantasia.wav');
L=length(y); %Longitud total del audio
while (index<L)
    if(index+0.025*Fs <L)
        next_index=index+0.025*Fs-1;
    else
        next_index=L;
    end
    audio=audioread('fantasia.wav',[index,next_index]); %Señal s
    autocorrelacion=xcorr(audio);
    [maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
    nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
    vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo a 19 muestras mas
    matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
    vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
    b=(matriz_ro^-1)*vector_ro; %Calculo los b
    g=sqrt(autocorrelacion(max_pos)-(transpose(b))*vector_ro); %Calculo la ganancia
    [h, w]=freqz(g,[1;((-1)*b)],1024);
    for j=1:1024
        Hmatrix(j,i)=abs(h(j)); %Armo la matriz H
        wmatrix(j,i)=w(j); %Armo la matriz de frecuencias  
    end
    t(i)=index; %Armo el vector de tiempos
    index=index+0.010*Fs; %Me corro 10 ms 
    i=i+1; %Me corro a la próxima columna de H
end
figure(5)
surface(t,wmatrix*Fs/(2*pi),log(Hmatrix))
colormap jet
colorbar
shading interp
view(2)
title('Superficie completa de las envolventes')
xlabel('Samples')
ylabel('Frecuencia')
figure(6)
colormap jet
spectrogram(y,1024,512,'yaxis')
title('Espectrograma de la señal original')


%%
%LPC APLICADO A CODIFICACION (item 5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ERROR CODIFICACION LPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,Fs] = audioread('fantasia.wav');
L=length(y); %Longitud total del audio
z0=zeros(20,1); %Condiciones iniciales para el filtro
audio_read= audioread('fantasia.wav'); %Leo el audio completo y lo completo con ceros
for m=L+1:L+365
    audio_read(m,1)=0;
end
i=1; %Posición en la columna de H
index=1; %Posición de lectura del audio
next_index=index+0.025*Fs-1; %Posicion del final de la ventana
longitud_ceros=length(audio_read);

while (next_index<=longitud_ceros)
    n=1;
    %Señal s=audio, tomo ventanas de 25 mseg separadas 10 mseg
    for m=index:next_index
        audio(n,1)=audio_read(m,1);
        n=n+1;
    end    
    autocorrelacion=xcorr(audio);
    [maximo,max_pos]=max(autocorrelacion); %Busco el valor donde la autocorrelacion es maxima
    nro=19; %Cantidad de muestras-1, ya que necesito 20 muestras
    vector_T=autocorrelacion([max_pos:max_pos+nro]); %Tomo del valor maximo a 19 muestras mas
    matriz_ro = toeplitz(vector_T); %Le paso la primer fila a la matriz de Toeplitz y la calculo
    vector_ro= autocorrelacion([max_pos+1 : max_pos+1+nro]); %Armo el vector ro
    b=[1 ;(-1)*(matriz_ro^-1)*vector_ro]; %Calculo los coeficientes b
    
    %Comienzo a filtrar 10 mseg de señal de la muestra tomada.
    audio_10mseg=audio(1:0.010*Fs); %Tomo 10 mseg de señal
    [error, z0]=filter(b,1,audio_10mseg,z0); %Calculo mi señal error
    longitud_error=length(error);
    longitud_b=length(b);
    
    %Armo dos matrices, la primera con los errores, y la segunda con los
    %coeficientes b
    for j=1:longitud_error
        error_matrix(j,i)=error(j);
    end
    
    for j=1:longitud_b
        b_matrix(j,i)=b(j);
    end
    
    t(i)=index; %Armo el vector de tiempos
    index=index+0.010*Fs; %Me corro 10 ms 
    next_index=index+0.025*Fs-1; %Tomo 25 ms de señal
    i=i+1; %Me corro a la próxima columna de error_matrix y b_matrix
end

%Empiezo a reconstruir la señal, primero sin redondear error y luego
%redondeando
longitud_10ms=0.01*Fs;
cond_inicio=zeros(20,1);
bits=4;

for h=1:4  
    if(h==1)
        for j=1:i-1
            [y_estimada,cond_inicio]=filter(1,b_matrix(:,j),error_matrix([1:longitud_10ms],j),cond_inicio);
            if j==1
                salida=y_estimada;
            else
                salida=[salida;y_estimada];
            end
        end
        %Tomo muestras del tamaño de la señal original
        for m=1:L
            audio_final(m,1)=salida(m,1);
        end
        error2=y-audio_final;
        estimaciones=audio_final;
    
    else
        %Matriz del error normalizada
        error_matrix_n=error_matrix;
        min_value=min(error_matrix_n(:));
        error_matrix_n=error_matrix_n-min_value;
        max_value=max(error_matrix_n(:));
        error_matrix_n=error_matrix_n/max_value;
        cond_inicio=zeros(20,1);       
        %Se redondea el error de prediccion a punto fijo de 4, 6 y 8 bits
        precision=2^bits-1;
        error_matrix_n=precision*error_matrix_n;
        error_matrix_n=round(error_matrix_n);
        %Asi se recibiria la señal, entonces, previo a recuperarla:
        error_matrix_n=error_matrix_n*max_value/precision;
        error_matrix_n=error_matrix_n+min_value; 
        
        %Se recupera la señal
            for j=1:i-1
                [y_estimada,cond_inicio]=filter(1,b_matrix(:,j),error_matrix_n([1:longitud_10ms],j),cond_inicio);
                if j==1
                    salida=y_estimada;
                else
                    salida=[salida;y_estimada];
                end
            end
       %Tomo muestras del tamaño de la señal original
       for m=1:L
        audio_final(m,1)=salida(m,1);
       end
       estimaciones=[estimaciones audio_final];
       error_nuevo=y-audio_final;
       error2=[error2 error_nuevo];
       bits=bits+2;
    end
end
%En la matriz error2 se graban como columna los errores, y en la matriz
%estimaciones se guardan como columna las señales estimadas en cada caso

%Plot de la señal reconstruida sin redondeo superpuesta con la original
reconstruccion_sinredondeo=estimaciones(:,1);
figure(7)
hold on
plot(y,'b')
plot(reconstruccion_sinredondeo,'r')
legend('Señal original','Señal estimada','Location','northeast')
title('Señales original y reconstruida superpuestas')
xlabel('Muestra')
ylabel('Amplitud')


%Plot de los errores
error_sinredondeo=error2([10000:10500],1);
error_4bits=error2([11000:11500],2);
error_6bits=error2([11000:11500],3);
error_8bits=error2([11000:11500],4);
figure(8)
hold on
plot(error_sinredondeo,'b')
plot(error_4bits,'r')
plot(error_6bits,'g')
plot(error_8bits,'k')
title('Error entre la señal real y su estimacion')
xlabel('Muestra')
ylabel('Amplitud')
legend('Error sin redondeo','Error 4 bits','Error 6 bits','Error 8 bits','Location','southwest')

%Plot de las señales reconstruidas
reconstruccion_sinredondeo=estimaciones([11000:11500],1);
reconstruccion_4bits=estimaciones([11000:11500],2);
reconstruccion_6bits=estimaciones([11000:11500],3);
reconstruccion_8bits=estimaciones([11000:11500],4);
figure(9)
hold on
plot(reconstruccion_sinredondeo,'b')
plot(reconstruccion_4bits,'r')
plot(reconstruccion_6bits,'g')
plot(reconstruccion_8bits,'k')
title('Señales reconstruidas')
xlabel('Muestra')
ylabel('Amplitud')
legend('Reconstrucción sin redondeo','Reconstrucción 4 bits','Reconstrucción 6 bits','Reconstrucción 8 bits','Location','southwest')

