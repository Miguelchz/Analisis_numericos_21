clear all
clc

while true
disp('ANALISIS NUMERICO 2021')
disp( '')
disp('unidad [1]')
disp('unidad [2]')
disp('unidad [3]')
disp('unidad [4]')
disp('unidad [5]')
disp('salir [otro numero]')
opcion=input('Ingrese una opcion:');
 
switch(opcion)
case 1
   
    %
clear all
clc

while true
disp('UNIDAD # 1'),
disp('SERIE DE TAYLOR [1]')
disp('salir [otro numero]')
opcion=input('ingrese opcion:');
 
switch(opcion)
case 1
    %INICIA UNIDAD 1****************************************************
    disp('saludo2'),
   % function y=taylor1(Y,n)
Y=input('ingrese la funcion: ');
n=input('ingrese el orden de la serie: ');
syms x

Y=eval(Y);
y=0;
x=0;
for l=0:n
y=y+eval(diff(Y,l))*(sym('x')^l)/factorial(l);
end
%end
    
  
    otherwise
        break
 end
 
end
    
    %TERMINA UNIDAD 1*****************************************************
    
    case 2
    
    %INICIA UNIDAD 2*****************************************************
clear all
clc
while true
 fprintf('\n')   
disp('UNIDAD # 2')
disp('Biseccion      [1]')
disp('falsa posicion [2]')
disp('punto fijo     [3]')
disp('newton         [4]')
disp('secante        [5]')
disp('Bairstow       [6]')
disp('Muller         [7]')
disp('salir          [otro numero]')
opcion=input('ingrese opcion:');
 
switch(opcion)
case 1
    disp('metodo de biseccion'),
    %function [c,k]= biseccion(y,a,b,e)
    clc
    syms x
    y=input('Ingrese la funcion: ')
    a=input('Ingrese limite inferior: ')
    b=input('Ingrese limite superior: ')
    e=input('Ingrese el error: ')
    fa=subs(y,x,a);
    fb=subs(y,x,b);
    k=0;
    if(fa*fb>0)
        disp('No esxites raiz!!!')
        %break
    end
    error=b-a;
    while error>e
    %while k<=6
        c=(a+b)/2;
        fc=subs(y,x,c);
        if(fc==0)
            a=c;
            b=c;
        elseif fa*fc<0
            b=c;
            fb=fc;
        else
            a=c;
            fa=fc;
        end
        error=b-a;
        k=k+1;
    end
    fprintf('El valor de la raiz es : %0.2f\n',c)
    fprintf('El numero de interacciones es : %d\n',k)

case 2
    disp('metodo de falsa posicion')
    clear,  
e = input('ingrese la funcion a evaluar: ');
f = inline (e)
 
x= -5:0.01:5;
y= f(x);
plot (x,y)
grid on
xa = input('ingrese el intervalo inferior: ');
xb = input('ingrese el intervalo superior: ');
t = input('ingrese el valor de la tolerancia: ');
vc = 0; ite =0.0000; error =100;
 
fprintf (' \t iteracion \t\txa \t\t vc \t\t xb \t error \n ')
 
while(error > t )
    
    vc=(xa*f(xb)-xb*f(xa))/(f(xb)-f(xa));
    disp([ite,xa,vc,xb,error])
    if(f(xa)*f(vc) < 0)
        xb=vc;
    else xa=vc;
    end
    error = abs(f(vc));
    ite= ite+1;
end
fprintf('la raiz se encontro en el punto:\n\t%f\n',vc)

case 3
    disp('metodo de punto fijo')
     %PuntoFijo(x0, es, imax, gx) %x0= valor inicial, es=porcentaje de error ,imax=iteraciones, gx=funcion despejada
     
     clear
     x0 = input('ingrese el valor inicial: ');
     es = input('ingrese el porcentaje de error: ');
     imax = input('cuantas iteraciones: ');
     gx = input('ingrese la funcion gx: ');
xr = x0;
iter = 0;
g = inline(gx);
do = 0;
while (do == 0)
    xrold = xr;
    xr = g(xrold);
    iter = iter + 1;
    if (xr ~= 0)
        ea = abs((xr - xrold)/xr)*100;
    end
    if ((ea < es) || (iter >= imax))
        break;
    end
end
disp('Resultado')
xr

case 4
    disp('metodo de newton raphson')
    clear
clc
disp('Metodo de newton Raphson')
syms x
f=input('introduce la funcion: ');
pi=input('introduzca el punto xi: ');
err=input('introduzca el porcentaje de error: ');
ezplot(f)
grid on
d=diff(f);% diff deriva una funcion
d=inline(d);
f=inline(f);
ea=100;%error aproximado
j=0;
while ea> err
    xi=pi-(f(pi)/d(pi));
    ea=abs(((xi-pi)/xi)*100);
    pi=xi;
    j=j+1;
end
fprintf('\nRaiz= %8.3f con %d Iteraciones ',pi,j);
 

case 5
    disp('metodo de la secante')
    clear
clc
disp('Metodo de la secante')
syms x
f=input('introduce la funcion: ','s');
x1=input('introduzca el punto xi-1: ');
x2=input('introduzca el punto xi: ');
err=input('introduzca el porcentaje de error : ');
ezplot(f)% dibuja la funcion
grid on
f=inline(f);
ea=100;
i=0;
fprintf('iteracion:           Raiz\n')
while ea>err
    xi=x2-(f(x2)*(x1-x2))/(f(x1)-f(x2));
    ea=abs((xi-x2)/xi*100);
    fprintf('%f           %8.4f\n',i,xi);
    x1=x2;
    x2=xi;
    i=i+1;
end
fprintf('\nRaiz de la funcion:%8.4f\ncalculada en %4fIteraciones\n',xi,i);

case 6
    disp('metodo de Bairstow')
    clear all; %limpio variables
clc; %Limpio pantalla
 
syms x % declaro variables simbolicas
 
t =[]; %declaro variables a usar
y =[]; %guarda primer division sintetica
e=[];  %permite reducir orden del polinomio
z=[];  %guarda segunda division sintetica
u=0;   %permite edentificar cada caso
n=2;   %permite asegurar el orden del polinomio se reduzca de dos en dos
g=0;
h=[];   %guarda salidas raices
o=[];   %guarda salidas raices
errora=1; %inicializo error
errorb=1;
tool=0; %inicializo tolerancia
f=input('\nIngrese la funcion: '); %pide ingresar la funcion
 
coeffs(f); %calcula coeficientes polinomio
t=sym2poly(f); %reorganiza polinomio
k=length(t); %calcula grado del polinomio
u=length(t); %variable por conveniencia, permite diferenciar si polinomio es mayor a grado 3 
p=length(t); % variable por conveniencia, se usa como contador 
x=k; % variable por conveniencia, permite imprimir cada caso
 
a=input('Ingrese punto a: '); %pide ingresar puntos (r,s) y tolerancia
b=input('Ingrese punto b: ');
tool=input('\nIngrese la tolerancia: ');
y(1)=1;% asegura que el primer termino de cada division sea 1
z(1)=1;
 
 
 
 
while t(1)~=1 %Este ciclo obliga a que el termino que acompa�a la variable de mayor grado deba ser igual a 1
    fprintf('\n\nerror el primer termino deber ser de valor 1, vuelva a ingresar polinomio...' )
    clear all; %limpio variables
    syms x % declaro variables simbolicas
    t =[]; %declaro variables a usar
    y =[];
    e=[];
    u=0;
    n=2;
    g=0;
    h=[];
    o=[];
    errora=1; %inicializo error
    errorb=1;
    tool=0; %inicializo tolerancia
    f=input('\nIngrese la funcion: '); %pide ingresar la funcion
    coeffs(f); %calcula coeficientes polinomio
    t=sym2poly(f); %reorganiza polinomio
    k=length(t); %calcula grado del polinomio
    u=length(t); %variable por conveniencia, permite diferenciar si polinomio es mayor a grado 3 
    p=length(t); % variable por conveniencia, se usa como contador 
    x=k; % variable por conveniencia, permite imprimir cada caso
    a=input('Ingrese punto a: '); %pide ingresar puntos y tolerancia
    b=input('Ingrese punto b: ');
    tool=input('\nIngrese la tolerancia: ');
    y(1)=1;
    z(1)=1;
end
 
 
if u==2 %si el polinomio a resolver es de grado 1 la raiz sale directamente 
    r1=-t(2);
    tool=1;
end
 
 
if  u==3 %Si el polinomio es de grado 2 se realiza el siguiente ciclo   
r1=(-t(2)+((t(2))^2-4*t(3))^(1/2))/2; %calcula raices
r2=(-t(2)-((t(2))^2-4*t(3))^(1/2))/2; 
end
 
 
while u>3 %si el polinomio ingresado es de grado 3 o mayor se realiza el siguiente ciclo
for i=1:k-1
    if i<2
    y(i+1)=a*t(i)+t(i+1);%permite calcular primer division sintetica
    else
    y(i+1)=a*y(i)+b*y(i-1)+t(i+1);
    end
end    
 
for i=1:k-2 %permite calcular segunda division sintetica
    if i<2
    z(i+1)=a*y(i)+y(i+1);
    else
    z(i+1)=a*z(i)+b*z(i-1)+y(i+1);
    end
end
y=fliplr(y);%reorganiza polinomios resultantes de las divisiones respectivas
z=fliplr(z);
 
da=(y(1)*z(3)-y(2)*z(2))/((z(2))^2-z(1)*z(3));%calculo pasos del metodo
db=(y(2)*z(1)-y(1)*z(2))/((z(2))^2-z(1)*z(3));
 
 
 
y=fliplr(y); %reorganiza polinomios resultantes de las divisiones respectivas
z=fliplr(z);
 
errora=da/a;%recalcula errores
errorb=db/b;
 
a=a+da; %calcula nuevos valores de a y b
b=b+db;
u=length(y);%calcula el grado del polinomio obtenido como cociente
 
 if tool>abs(errora) | tool>abs(errorb) & u==1 %si polinomio del cociente es de grado 1 
         r5=(-a/b);        
     end
 if tool>abs(errora) | tool>abs(errorb) & u==2 %si polinomio del cociente es de grado 2 
         r3=(a+(a^2+4*b)^(1/2))/2; 
         r4=(a-(a^2+4*b)^(1/2))/2;
     end
    if tool>abs(errora) | tool>abs(errorb) & u>=3 %si polinomio del cociente es de grado 3 o mayor          
         r1=(a+(a^2+4*b)^(1/2))/2; 
         r2=(a-(a^2+4*b)^(1/2))/2; 
         g=g+1;
         o(g)=r1; %guarda las raices 
         h(g)=r2; %guarda las raices
         for j=1:p-n %permite reducir orden del polinomio resultante
             e(j)=y(j);
         end
         y=e; %se asignan nuevos valores a variables necesarias
         z=[];
         z(1)=1;
         t=e;             
         k=length(y);
         u=length(y);
         n=n+2; % asegura que el grado del polinomio se disminuya de 2 en 2
         errora=1; %reinicia errores
         errorb=1;
         e=[];      
         end
end
 
fprintf('Las raices son: \n\n');
 
if x==2
 r1 
end
 
if x==3
 r1 
 r2
end
 
if x>3 & u==2
    o
    h
    r5
end
 
if x>3 & u==3
    o
    h
    r3
    r4
end
 
fprintf('\nLos valores de r y s son: %f',a);
fprintf('\n%f',b);

case 7
    disp('metodo de muller')
    close all; % cierra ventana
clear; % limpia variables
clc % limpia
%%
syms x
xr = input('xr = ');
h = input('h = ');
eps = input('eps = ');
maxit = input('maxit =');
fx = input('f(x) = '); % x^2 + 3*x + 2
x2 = xr;
x1 = xr + h*xr;
x0 = xr - h*xr;
f = inline(fx); % subs(fx,x,x0)
%%
iter = 0;
while(iter < maxit)
   iter = iter + 1; 
   h0 = x1 - x0;
   h1 = x2 - x1;
   d0 = (f(x1) - f(x0))/h0;
   d1 = (f(x2) - f(x1))/h1;
   a = (d1 - d0)/(h1 + h0);
   b = a*h1 + d1;
   c = f(x2);
   rad = sqrt(b*b - 4*a*c);
   if(abs(b+rad) > abs(b-rad))
       den = b +rad;
   else
       den = b - rad;
   end
   dxr = -2*c/den;
   xr = x2 + dxr;
   if(abs(dxr) < eps*xr)
       break;
   end
   x0 = x1;
   x1 = x2;
   x2 = xr;  
end
disp('xr')
disp(xr)
    otherwise
        clc
        disp('a salido de unidad 2')
        break
 end

end
    
    
    
    
    %TERMINA UNIDAD 2******************************************************
    case 3
    disp('UNIDAD #3')
    
     %INICIA UNIDAD 3******************************************************
     
     clear all
     clc
while true
    fprintf('\n')
       fprintf('UNIDAD # 3')
       fprintf('\n')
disp('interpolacion de la grange [1]')
disp('interpolacion de newton    [2]')
disp('diferencias Divididas      [3]')
disp('polinomio de hermite       [4]')
disp('salir [otro numero]')
opcion=input('ingrese opcion:');
 
switch(opcion)
case 1
    disp('interpolacion de la grange'),
    %function p=interpolacionlagrange(x,f,a)% vectores x y vectores f y el punto a
x=input('ingrese los valore de "x" :');
f=input('ingrese los valore de "y" :');
a=input('ingrese el punto "a" para interpolar :');
n=length (x);
syms t;
p=0;
plot (x,f,'*')
grid on;
for i=1:n
    L=1;
    for j=1:n
        if i~=j
            L=L*(t-x(j))/(x(i)-x(j));
        end
    end
    p=p+L*f(i);
end
p=expand(p)
hold on;
ezplot(p,[0,10])

t=a;
p=eval(p)
    case 2
    fprintf('interpolacion de newton')
    fprintf('\n')
   % function interpolaNewton
    %Interpolacion de newton
    clear;clc;
    %disp('metodos numericos');
    %disp('interpolacion');
    %disp('interpolacion');
    n=input('ingrese el grado del polinomio, n=');
    fprintf('Se necesitan %.0f puntos\n',n+1);
    disp('ingrese los puntos');
    for i=1:n+1
        fprintf('x%.0f=',i-1);
        X(i)=input(' ');
        fprintf('y%.0f=',i-1);
        Y(i)=input(' ');
    end
    DD=zeros(n+1);
    DD(:,1)=Y;
    for k=2:n+1
        for J=k:n+1
            DD(J,k)=[DD(J,k-1)-DD(J-1,k-1)]/[X(J)-X(J-k+1)];
        end
    end
    disp('La matriz de diferencias divididas es:');
    disp(DD);
    disp('El polinomio de newton es');
    syms x;
    polnew=DD(1,1);
    P=1;
    for i=1:n
        P=P*(x-X(i));
        polnew=polnew+P*DD(i+1,i+1);
    end
    polnew=expand(polnew);
    pretty(polnew);
    x=input('ingrese el valor de x a interpolar,x=');
    vi=eval(polnew);
    fprintf('el valor interpolado es %.2f\n',vi);
    hold on;
    ezplot(polnew,[X(1) X(n+1)]);
    plot(x,vi,'r+');
    case 3
    disp('diferencias Divididas')
    x=input('introduzca el vector en "x"');
y=input('introduzca el vector en "y"');
xi=input('introduzca el valor "xi" a interpolar :');
n=length(x);
b=zeros(n);
b(:,1)=y(:);

for j=2:n
    for i=1:n-j+1
        b(i,j)=(b(i+1,j-1)-b(i,j-1))/(x(i+j-1)-x(i))
    end
end

x1=1;
yi=b(1,1);
for j=1:n-1
    x1=x1.*(xi-x(j));
    yi=yi+b(1,j+i)*x1;
end
disp('el valor de "yi" interpolado es igual a'),disp (yi);

p=num2str(b(1,1));
xx=x*-1;
for j=2:n
    signo='';
    if b(1,j)>=0
        signo='+';
    end
    x1='';
    for i=1:j-1
        signo2='';
        if xx(i)>=0
            signo2='+';
        end
        x1=strcat(x1,'*(x)',signo2,num2str(xx(i)),')');
    end
    p=strcat(p,signo,num2str(b(1,j)),x1);
end
disp('el polinomio de interpolacion de newton es :');
p=inline(p)
plot(x,y);
grid on 


    case 4
    disp('polinomio de hermite')
    X=input('Ingrese los valores de x='); % en forma de vector
Y=input('Ingrese los valores de f(x)='); % en forma de vector
DF=input('Ingrese los valores de la derivada de f(x)='); % en forma de vector
x=input('Ingrese el valor a interpolar ='  );
n=length(X);
Q=zeros(2,n);
for i=1:n
z(2*i-1)=X(i);
z(2*i)=X(i);
Q(2*i-1,1)=Y(i);
Q(2*i,1)=Y(i);
Q(2*i,2)=DF(i);
if i~=1
Q(2*i-1,2)=(Q(2*i-1,1)-Q(2*i-2,1))/(z(2*i-1)-z(2*i-2));
end
end
for i=3:2*n
for j=3:i
Q(i,j)=(Q(i,j-1)-Q(i-1,j-1))/(z(i)-z(i-j+1));
end
end
syms x
Fx=Q(1,1);
%Diferencias divididas
for p=1:numel(X)-1
L=1;
%Multiplicaci�n de los polinomios
for k=1:p
L=L*(x-X(k));
end
Fx=Fx+L*Q(p+1,p+1);
end
%Aproximacion del Polinomio resultante
val=eval(Fx);
disp(val);
    otherwise
        clc
        disp('a salido de la unidad 3')
        break
 end
 
end

    
    %TERMINA UNIDAD 3******************************************************
    
    case 4
    disp('UNIDAD #4')
    disp(' ')
     %INICIA UNIDAD 4******************************************************
     %menu para las opciones de la unidad 4
     while true
     disp('Derivaci�n                    [1]')
     disp('Integraci�n                   [2]')
     disp('Extrapolaci�n de Richardson   [3]')
     disp('Integraci�n por Rosemberg     [4]')
     disp('Regresar                      [5]')
     disp(' ')
     opcion=input('Ingrese una opcion:');
     
     %switch para las opciones
     %---------------------------------------------------------------------
     switch(opcion)
     %Desde aqui van todos los case
         case 1
             %empieza el codigo para la opcion de derivacion
clear all
fprintf('\n');
disp('----Derivaci�n num�rica----');
%obteniendo los valores
%para poder usar la x de las funciones
syms x
%obteniendo la funcion que se desea derivar
g=input('Ingrese la funci�n: ');
%obteniendo el valor de h
hIng=input('Ingrese el valor de h: ');
%obteniendo el valor de que se desea evaluar
ax=input('Ingrese el valor a evaluar: ');
format long
f=inline(g,'x') 
%Inicio del while para todo el menu
while true
 fprintf('\n')   
disp('Opciones disponibles: ');
disp(' ');
disp('Diferencias Finitas           [1]')
disp('F�rmula de los 3 puntos       [2]')
disp('F�rmula de los 5 puntos       [3]')
disp('Derivadas de orden superior   [4]')
disp('Salir                         [5]')
disp(' ')
opcion=input('Ingrese una opcion:');

%empiezan opciones para cada seleccion
switch(opcion)
%opciones para diferencias finitas    
case 1
        fprintf('\n');
        disp('----Diferencias finitas----');
        while true
            disp('Hacia adelante    [1]')
            disp('Hacia atras       [2]')
            disp('Centrada          [3]')
            disp('Regresar          [4]')
            disp(' ')
            opcion=input('Ingrese una opcion:');
            switch(opcion)
                %hacia adelante
                case 1
                    while true
                       disp('Primera diferencia    [1]')
                       disp('Segunda diferencia    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %hacia adelante primera diferencia    
                                    r=(f(ax+hIng)-f(ax))/(hIng);
                                    fprintf('El resultado hacia adelante primera diferencia es: \n\t%f\n',r)
                                    disp(' ');
                                    break                                    
                                    case 2   
                                    %hacia adelante segunda diferencia
                                    r=(-f(ax+2*hIng)+4*f(ax+hIng)-3*f(ax))/(2*hIng);
                                    fprintf('El resultado hacia adelante segunda diferencia es: \n\t%f\n',r)
                                    disp(' ');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break                                   
                                    end
                    end                 
                    %hacia atras
                case 2
                     while true
                       disp('Primera diferencia    [1]')
                       disp('Segunda diferencia    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %hacia atras primera diferencia    
                                    r=(f(ax)-f(ax-hIng))/hIng;
                                    fprintf('El resultado hacia atras primera diferencia es: \n\t%f\n',r)
                                    disp(' ');
                                    break         
                                    case 2   
                                    %hacia atras segunda diferencia
                                    r=(3*f(ax)-4*f(ax-hIng)+f(ax-2*hIng))/(2*hIng);
                                    fprintf('El resultado hacia atras segunda diferencia es: \n\t%f\n',r)
                                    disp(' ');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break
                                    end
                     end
                %centrada
                case 3
                     while true
                       disp('Orden 2    [1]')
                       disp('Orden 4    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %centrada primera diferencia    
                                    r=(f(ax+hIng)-f(ax-hIng))/(2*hIng);
                                    fprintf('El resultado de la centrada orden 2 es: \n\t%f\n',r)
                                    disp(' ');
                                    break                    
                                    case 2   
                                    %centrada segunda diferencia
                                    r=(-f(ax+2*hIng)+8*f(ax+hIng)-8*f(ax-hIng)+f(ax-2*hIng))/(12*hIng);
                                    fprintf('El resultado de la centrada orden 4 es:\n \n\t%f\n',r)
                                    disp(' ');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break                                    
                                   end
                     end 
                case 4
               %regesa al menu
                break
            end
        %aqui termina el case para diferencias finitas    
        end
%aqui comienza la formula de los 3 puntos
case 2
        fprintf('\n');
        disp('----F�rmula de los 3 puntos----');
        while true
            disp('Formula 1         [1]')
            disp('Formula II        [2]')
            disp('Regresar          [3]')            
            disp(' ')
            opcion=input('Ingrese una opcion:');
            switch(opcion)
                %Formula I
                case 1
                     r=(f(ax+hIng)-f(ax-hIng))/(2*hIng);
                     fprintf('El resultado de la formula I es: \n\t%f\n',r)
                     disp(' ');
                  break                  
                    %Formula II
                case 2
                     r=(((-3*f(ax))+4*f(ax+hIng)-f(ax+2*hIng))/(2*hIng));
                     fprintf('El resultado de la formula II es: \n\t%f\n',r)
                     disp(' ');
                    break
                case 3
               %regesa al menu
                break
            end
        %aqui termina el case para la formula de los 3 puntos
        end
%aqui comienzan las opciones de la formula de los 5 puntos
case 3
        fprintf('\n');
        disp('----F�rmula de los 5 puntos----');
        while true
            disp('F�rmula I      [1]')
            disp('F�rmula II     [2]')
            disp('F�rmula III    [3]')
            disp('F�rmula VI     [4]')
            disp('F�rmula III    [5]')
            disp('Regresar       [6]')
            disp(' ')
            opcion=input('Ingrese una opcion:');
            switch(opcion)
                %F�rmula 1
                case 1
                     r=((1)/(12*hIng))*((-25*f(ax))+(48*(f(ax+hIng)))-(36*(f(ax+2*hIng)))+(16*(f(ax+3*hIng)))-(3*(f(ax+4*hIng))));
                     fprintf('El resultado de la formula I es: \n\t%f\n',r)
                     disp(' ');
                   break
                    %F�rmula II
                case 2
                     r=((1)/(12*hIng))*((-3*f(ax-hIng))-(10*(f(ax)))+(18*(f(ax+hIng)))-(6*(f(ax+2*hIng)))+(f(ax+3*hIng)));
                     fprintf('El resultado de la formula II es: \n\t%f\n',r)
                     disp(' ');
                    break
                    %Formula III
                case 3
                     r=((1)/(12*hIng))*((f(ax-2*hIng))-(8*(f(ax-hIng)))+(8*(f(ax+hIng)))-((f(ax+2*hIng))));
                     fprintf('El resultado de la formula III es: \n\t%f\n',r)
                     disp(' ');
                    break
                    %IV 
                case 8
                     r=((1)/(12*hIng))*((4*f(ax-3*hIng))+(6*(f(ax+2*hIng)))-(8*(f(ax-hIng)))+(34*(f(ax)))+(3*(f(ax+hIng)))+(34*(f(ax+2*hIng))));
                     fprintf('El resultado de la formula IV es: \n\t%f\n',r)
                     disp(' ');
                    break
                %formula v
                case 5
                     r=((1)/(12*hIng))*((f(ax-4*hIng))-(3*(f(ax-3*hIng)))+(4*(f(ax-2*hIng)))-(36*(f(ax-hIng)))+(25*(f(ax))));
                     fprintf('El resultado de la formula V es: \n\t%f\n',r)
                     
                     disp(' ');
                    break
                case 6
                    %regresa al menu
                    break
            end
        %aqui termina el case para la formulas de 5 puntos
        end   
%empizan las opciones para derivadas de orden superior        
case 4
        fprintf('\n');
        disp('----Derivadas de orden superior----');
        while true
            disp('Hacia adelante    [1]')
            disp('Hacia atras       [2]')
            disp('Centrada          [3]')
            disp('Regresar          [4]')
            disp(' ')
            opcion=input('Ingrese una opcion:');
            switch(opcion)
                %hacia adelante
                case 1
                    while true
                       disp('Primera diferencia    [1]')
                       disp('Segunda diferencia    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %hacia adelante primera diferencia    
                                    r=(f(ax+2*hIng)-2*f(ax+hIng)+f(ax))/(hIng^2);
                                    r1=(f(ax)-3*f(ax-hIng)+3*f(ax-2*hIng)-f(ax-3*hIng))/(hIng^3);
                                    r2=(f(ax)-4*f(ax-hIng)+6*f(ax-2*hIng)-4*f(ax-3*hIng)+f(ax-4*hIng))/(hIng^4);
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break                                   
                                    case 2   
                                    %hacia adelante segunda diferencia
                                    rr=(-f(ax+3*hIng)+4*f(ax+2*hIng)-5*f(ax+hIng)+2*f(ax))/(hIng^2);
                                    r1=(-3*f(ax+4*hIng)+14*f(ax+3*hIng)-24*f(ax+2*hIng)+18*f(ax+hIng)-5*f(ax))/(2*(hIng^3));
                                    r2=(-2*f(ax+5*hIng)+11*f(ax+4*hIng)-24*f(ax+3*hIng)+26*f(ax+hIng)+3*f(ax))/(hIng^4);
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break                                 
                                    end
                    end                     
                    %hacia atras
                case 2
                     while true
                       disp('Primera diferencia    [1]')
                       disp('Segunda diferencia    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %hacia atras primera diferencia    
                                    r=(f(ax)-2*f(ax-hIng)+f(ax-2*hIng))/(hIng^2);
                                    r1=(f(ax)-3*f(ax-hIng)+3*f(ax-2*hIng)-f(ax-3*hIng))/(hIng^3);
                                    r2=(f(ax)-4*f(ax-hIng)+6*f(ax-2*hIng)-4*f(ax-3*hIng)+f(ax-4*hIng))/(hIng^4);
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break                                    
                                    case 2   
                                    %hacia atras segunda diferencia
                                    r=(2*f(ax)-5*f(ax-hIng)+4*f(ax-2*hIng)-f(ax-3*hIng))/(hIng^2);
                                    r1=(5*f(ax)-18*f(ax-hIng)+24*f(ax-2*hIng)-14*f(ax-3*hIng)+3*f(ax-4*hIng))/(hIng^3);
                                    r2=(3*f(ax)-14*f(ax-hIng)+26*f(ax-2*hIng)-24*f(ax-3*hIng)+11*f(ax-4*hIng)-2*f(ax-5*hIng))/(hIng^4);
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break                                   
                                    end
                     end                   
                    %centrada
                case 3
                     while true
                       disp('Primera diferencia    [1]')
                       disp('Segunda diferencia    [2]')
                       disp('Regresar              [3]')
                       disp(' ')
                                opdif=input('Ingrese una opcion:');
                                    switch(opdif)
                                    case 1
                                    %centrada primera diferencia    
                                    r=(f(ax+hIng)-2*f(ax)+f(ax-hIng))/(hIng^2);
                                    r1=(f(ax+2*hIng)-2*f(ax+hIng)+2*f(ax-hIng)-f(ax-2*hIng))/((2)*(hIng^3));
                                    r2=(f(ax+2*hIng)-4*f(ax+hIng)+6*f(ax)-4*f(ax-hIng)+f(ax-2*hIng))/(hIng^4);
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break                                 
                                    case 2   
                                    %centrada segunda diferencia
                                    r=(-f(ax+2*hIng)+16*f(ax+hIng)-30*f(ax)+16*f(ax-hIng)-f(ax-2*hIng))/((12)*(hIng^2));
                                    r1=(-f(ax+3*hIng)+8*f(ax+2*hIng)-12*f(ax+hIng)+12*f(ax-hIng)-8*f(ax-2*hIng)+f(ax-3*hIng))/((8)*(hIng^3));
                                    r2=(-f(ax+3*hIng)+12*f(ax+2*hIng)-39*f(ax+hIng)+56*f(ax)-39*f(ax-hIng)+12*f(ax-2*hIng)-f(ax-3*hIng))/((6)*(hIng^4));
                                    fprintf('El resultado es:\n');
                                    disp(' ');
                                    fprintf('f''''(x): ');
                                    fprintf('%f',r);
                                    fprintf('\nf''''''(x): ');
                                    fprintf('%f',r1);
                                    fprintf('\nf''''''''(x): ');
                                    fprintf('%f',r2);
                                    fprintf('\n------------------\n');
                                    break
                                        case 3 
                                        %regresa al menu
                                            break                                    
                                   end
                     end                   
                case 4
               %regesa al menu
                break
            end
        %aqui termina el case para derivadas de orden superior
        end
        case 5
    break
end
        %aqui termina el case para diferencias finitas     end
%aqui terminan las opciones para todos los metodos
end
%fin del while del menu

             
             %termina el codigo para la derivacion
             break
         case 2 
              %empieza el codigo para la opcion de integracion
clear all
fprintf('\n');
disp('----Integraci�n num�rica----');
%obteniendo los valores
syms x
%funcion
g=input('Ingrese la funci�n: ');
%para poder usar la x de las funciones
%limite inferior
a=input('Ingrese el limite inferior: ');
%limite superior
b=input('Ingrese el limite superior: ');
f=inline(g,'x');
%Inicio del while para todo el menu
while true
fprintf('\n')   
disp('Opciones disponibles: ');
disp(' ');
disp('Regla del trapecio   [1]')
disp('Regla simpson 1/3    [2]')
disp('Regla simpson 3/8    [3]')
disp('Salir                [4]')
disp(' ')
opcion=input('Ingrese una opcion:');
switch(opcion)
case 1
%regla del trapecio  
r=(b-a)*((f(a)+f(b))/2);
r=(b-a)*((f(a)+f(b))/2);
fprintf('El resultado es :%f',r)
disp(' ');
 
  
%finaliza la regla del trapecio    
break
case 2
    %regla simpson 1/3
h=(b-a)/2;
r=(b-a)*((f(a)+4*f(h)+f(b))/6);
fprintf('El resultado es :%f',r)
    %finaliza simpson 1/3
break
case 3 
        %regla simpson 3/8
h=(b-a)/3;
ha=zeros(1,3);
ha(1)=a;
for k=2:3
  ha(k)=ha(k-1)+h;
end
r=(b-a)*((f(a)+3*f(ha(2))+3*f(ha(3))+f(b))/8);
fprintf('El resultado es :%f',r)
        %finaliza regla de simpson 3/8
break
case 4 
        %regresar
break
        
 %end del switch principal       
end       
%finaliza el while de todo el menu    
end
             
             %termina el codigo para la integracion
             break
         case 3 
             %empieza el codigo para la extrapolacion de richardson
             
             %termina el codigo para la extrapolacion de richardson
             break
         case 4 
             %empieza el codigo para la integracion de rosemberg
clear all
fprintf('\n');
disp('----Integracion de Romberg----');
%obteniendo los valores
%para poder usar la x de las funciones
syms x
%obteniendo la funcion que se desea derivar
f=input('Ingrese la funci�n: ');
a = input('Ingrese el limite inferior, a:  ');
b = input('Ingrese el limite superior, b:  ');
%obteniendo el valor de h
n=input('Ingrese el nivel ');
%obteniendo el valor de que se desea evaluar

format long
f=inline(f,'x') 
 
 h = b-a;
 r = zeros(2,n+1);
 r(1,1) = (f(a)+f(b))/2*h;
 fprintf('\nResultados de la integracion de Romberg:\n');
 fprintf('\n %11.8f\n\n', r(1,1));

 for i = 2:n
    sum = 0;
    for k = 1:2^(i-2)
       sum = sum+f(a+(k-0.5)*h);
    end
    r(2,1) = (r(1,1)+h*sum)/2;
   
    for j = 2:i
       l = 2^(2*(j-1));
       r(2,j) = r(2,j-1)+(r(2,j-1)-r(1,j-1))/(l-1);
    end

    for k = 1:i
       fprintf(' %11.8f',r(2,k));
    end
  
    fprintf('\n\n');
    h = h/2;
    for j = 1:i
       r(1,j) = r(2,j);
    end
 end
             %termina el codigo para la integracion de rosemberg
             break
         case 5 
            %regresar al menu
             break
     %Hasta aqui deben de llegar todos los case
     %termina el switch de las opciones     
     end    
     %---------------------------------------------------------------------
     
     %terminan las opciones de la unidad 4
     end
    %TERMINA UNIDAD 4******************************************************
    
    case 5
    disp('UNIDAD #5')
     %INICIA UNIDAD 5******************************************************
%inician las opcines para la unidad 5    
while true
disp('Euler                    [1]')
disp('Runge Kutta              [2]')
disp('Multipasos               [3]')
disp('Salir                    [4]')
opcion=input('Ingrese una opcion:');
switch(opcion)
case 1
    %Euler
    syms x
syms y %variable simbolica
f=inline(input('ingrese la derivada: ','s'));
x=input('ingrese el valor de x: '); %valor inicial
xf=input('ingrese el valor final de x: ');%valor ultimo
y=input('ingrese el valor de y: ');% el valor de variable y
h=input('ingrese el valor de h: '); %tamanio de paso
n=(xf-x)/h
clc
disp('x(n) y(n)');

for i=1:n+1
    y1= feval(f,x,y);
    hy1=h*y1;
    fprintf('\n%0.1f %0.4f',x,y);
    y=y+hy1;
    x=x+h;
end
break
case 2
    %Runge Kutta
    clear
clc
format long

F=input('ingrese la funcion del lado derecho f(x,y)= ','s');
funcion=inline(F, 'x','y');%crea la funcion en matlab
a=input('ingrese el valor de a: ');
b=input('ingrese el valor de b: ');
y0=input('ingrese el valor de y0: ');
h=input('ingrese el valor de h :');%tamanio de paso

N=(b-a)/h;%numero de puntos 

%matricez
x=zeros(1,N+1); %vector filas
y=zeros(1,N+1);
k1=zeros(1,N);
k2=zeros(1,N);
k3=zeros(1,N);
k4=zeros(1,N);
y(1)=y0; % define primer elemento de yi

%iteraciones
for i=1:N+1
    x(i)=a+(i-1)*h;
end
%genera los ki del metodo
for j=1:N
    k1(j)=funcion(x(j),y(j));
    k2(j)=funcion(x(j)+h/2,y(j)+h*k1(j)/2);
    k3(j)=funcion(x(j)+h/2,y(j)+h*k2(j)/2);
    k4(j)=funcion(x(j)+h,y(j)+h*k3(j));
    y(j+1)=y(j)+(h/6)*(k1(j)+2*k2(j)+2*k3(j)+k4(j));
end
%muestra los valores
fprintf('     i           xi               yi \n')
for k=1:N+1
    fprintf('%6.0f %12.12f %12.12f\n',k,x(k),y(k))
end
%grafica

plot(x,y,'ro','Linewidth',2)
title('grafica de iteraicones')
xlabel('eje x')
ylabel('eje y')
legend('(xi,yi)')
break
case 3 
    %multipasos
%f=input('introduzca la funcion derivada: ');
f=@(t,y) (y-t^2+1);%cambiar la funcion (y-t^2+1) por caulquier otra

a=input('introduzca el punto final izquierdo, a: ');
b=input('introduzca el punto final derecho, b: ');
n=input('introduzca el # de sub-intervalos, n: ');
alpha=input('introduzca la condicion inical: ');

h=(b-a)/n;
t(1)=a;
w(1)=alpha;
fprintf('   t          w\n');
fprintf('%5.4f  %11.8f\n', t(1), w(1));

for i=1:3
    t(i+1)= t(i)+h;
    k1 = h*f(t(i),w(i));
    k2 = h*f(t(i)+0.5*h,w(i)+0.5*k1);
    k3 = h*f(t(i)+0.5*h,w(i)+0.5*k2);
    k4 = h*f(t(i+1),w(i)+k3);
    w(i+1) = w(i)+(k1+2.0*(k2+k3)+k4)/6.0;
    fprintf('%5.4f  %11.8f\n',  t(i+1), w(i+1));
end

for i= 4:n
    t0 = a+i*h;
    part1=55.0*f(t(4),w(4))-59.0*f(t(3),w(3))+37.0*f(t(2),w(2));
    part2 = -9.0*f(t(1),w(1));
    w0 = w(4)+h*(part1+part2)/24.0;
    part1 = 9.0*f(t0,w0)+19.0*f(t(4),w(4))-5.0*f(t(3),w(3))+f(t(2),w(2));
    w0=w(4)+h*(part1)/24.0;
    fprintf('%5.4f  %11.8f\n',t0,  w0);
    for j=1:3
        t(j)= t(j+1);
        w(j)= w(j+1);
    end
    t(4) = t0;
    w(4)= w0;
end
break
case 4 
    %salir
break
   
end    

%terminan las opcines para la unidad 5
end
    
    %TERMINA UNIDAD 5******************************************************
    
    
    otherwise
        dispse
        clc('a salido de los metodos numericos!')
        break
 end
 
end
