function MDFOnda(Uzero,Ul,Uini,dUini,c,l,tmax,n,m,k)
%clear all
%close all
clc

syms x t1 n1

%estavel qdo h/k <= 1/c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            valores para teste

%Formula do MDF para eq da onda
%u(i,j+1)=mu^2*u(i+1,j)+2*(1-mu^2)*u(i,j)+mu^2*u(i-1,j)-u(i,j-1)

% l=1;              %Tamanho da corda
% tmax=0.5;            %tempo maximo de analise
% n=5;               %qtd de nó em x
% m=10;               %qtd de nó em t
% c=1;                %%constante da eq da onda
%Uzero=0              %condição de contorno em x=0
%Ul=0                 %condição de contorno em x=l
%Uini = sin(pi*x)     %condição inicial
%g = 0                %condição inicial derivada de u
%k = 200              %numero da soma parcial da serie de fourier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Solução analitica via serie de Fouriercn=(2/l)*int(f*sin((n1*pi*x)/l),x,0,l);
cn=(2/l)*int(Uini*sin((n1*pi*x)/l),x,0,l);
dn=(l/(n1*c*pi))*(2/l)*int(dUini*sin((n1*pi*x)/l),x,0,l);


a1=cn*sin((n1*pi*x)/l)*cos((n1*c*pi*t1)/l);
a2=dn*sin((n1*pi*x)/l)*sin((n1*c*pi*t1)/l);



solucao=((symsum(a1,n1,1,k))+(symsum(a2,n1,1,k)));


k=l/n;         %% passo em x
h=tmax/m;      %%passo em t

mu=(h*c)/k;
mu

x1=(0:k:l);
t=(0:h:tmax);


a=length(x1); %% quantidade de nós de x
b=length(t);  %% quantidade de instantes de tempo

U=zeros(a,b); %%Matriz onda cada coluna representa um instante de tempo e cada linha representa um nó em x

q=1;

% Calculo da primeira coluna da matriz U, que se refere a condição inicial
% da corda no tempo t=0
for i=0:k:l

    U(q,1)=double((subs(Uini,x,i)));
    
    q=q+1;
end

%Condição de Contorno
U(1,:)=Uzero;
U(a,:)=Ul;

q=2;

% Calculo da segunda coluna da matriz U

for i=k:k:(l-k)
    U(q,2)=((mu^2)/2)*(U(q+1,1)+U(q-1,1))+(1-mu^2)*(U(q,1))+h*(double((subs(dUini,x,i))));
    q=q+1;
end



 %calculo das outras colunas de U onde a coluna j é o indice do tempo e i o
 %indice do x
for j=2:b-1
   
   for i=2:(a-1)
            U(i,j+1)=mu^2*U(i+1,j)+2*(1-mu^2)*U(i,j)+mu^2*U(i-1,j)-U(i,j-1);
    end
   
end


count=1;
%% S é um vetor que tem os valores da solução analitica para cada nó de x
for i=0:h:(tmax)
    
    aux=subs(solucao,t1,i);
    q=1;
    for j=0:k:l
        S(count,q)=double(subs(aux,x,j))  ;
        q=q+1;
    end
    count=count+1;
   
end



U=U'
figure(1);
%%Definindo os limites maximos e minimo de altura da figura
ymax=max(max(U))+3;
ymin=min(min(U))-3;

%hold on
for i=1:b
    
    plot(x1,U(i,:),x1,S(i,:),'--','LineWidth',2)

    axis([0,l ,ymin, ymax])
    legend('Diferenças Finitas','Serie de Fourier')
    title('Comparação das Solução');
    ylabel('Deslocamento');
    xlabel('x');
    M(i) = getframe;
    
end  

    