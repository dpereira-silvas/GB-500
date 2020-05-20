function mdfCalor(Uini,Uzero,Ul,alfa,l,dx,dt,tmax,k)

clc
%Formula do metodo da diferen�as Finitas para eq do calor 1d
%u(i,j+1)=u(i,j)+lambda(u(i+1,j)-2*u(i,j)+u(i-1,j))
%lambda = (alfa^2)*(dt/(dx^2))
%estav�l Quando lambda <= 1/2


syms n1 x t1 %cria��o de variaveis simbolicas

%Dados de Entrada

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               valores para teste

% alfa = 1; Constante da eq do calor
% dx=5; espa�amento dos n�s em x
% dt=5; espa�amento dos instantes de tempo
% l=50; comprimento da barra
% tmax=101; tempo maximo de analise
% Uini=60-2*x; condi��o inicial
% Uini=20;
% Uzero=0 condi��o de contorno em x=0
% Ul=0  condi��o de contorno em x=l
%k = 200              %numero da soma parcial da serie de fourier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Solu��o analitica via serie de Fourier

cn=(2/l)*int((Uini-(Uzero+((Ul-Uzero)/l)*x))*sin((n1*pi*x)/l),x,0,l);

f=cn*sin((n1*pi*x)/l)*exp((-(n1^2+pi^2+alfa^2)/l^2)*t1);

%soma Parcial somatorio em k
% para k=200
solAnal=(Uzero+((Ul-Uzero)/l)*x)+(symsum(f,n1,1,k));



lambda = (alfa^2)*(dt/(dx^2));

x1=(0:dx:l);
t=(0:dt:tmax);


m=length(x1); %% quantidade de n�s de x
k=length(t); %% quantidade de instantes de tempo

U=zeros(m,k);%%Matriz onda cada coluna representa um instante de tempo e cada linha representa um n� em x



%condi��o inicial
q=1;
% Calculo da primeira linha da matriz U, que se refere a condi��o inicial
% da barra no tempo t=0
for i=0:dx:l

    %U(i,1)=20;
    U(q,1)=double((subs(Uini,x,i)));
    q=q+1;
end

%Condi��o de Contorno
 U(1,:)=Uzero; %U(x=0,t)
 U(m,:)=Ul;    %U(x=l,t)

 %calculo das outras linhas de U onde a coluna j � o indice do tempo e i o
 %indice do x
for j=2:k
   
   for i=2:(m-1)
            U(i,j)= ( U(i,j-1) + ( lambda*( U( i+1,j-1 )- ( 2*U(i,j-1) ) + U(i-1,j-1) ) ) );
    end
   
end

a=U';
count=1;
%% S � um vetor que tem os valores da solu��o analitica para cada n� de x
for i=0:dt:tmax
    
    b=subs(solAnal,t1,i);
    q=1;
    for j=0:dx:l
     S(count,q)=double(subs(b,x,j));
         q=q+1;
    end
    count=count+1;
end





figure(1);
%%Definindo os limites maximos e minimo de altura da figura
ymax=max(max(U))+3;
ymin=min(min(U))-3;

%hold on

for i=1:k
        plot(x1,a(i,:),x1,S(i,:),'--','LineWidth',2);
        axis([0,l ,ymin, ymax])
        text(5,15,['t= ',num2str(i)]);
        legend('Diferen�as Finitas','Serie de Fourier')
        ylabel('Temperatura')
        xlabel('x')
        M(i) = getframe;
     
end
