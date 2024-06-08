function [U] = PRNGNormal(n)
%PRNGUNIFORMCONG Questa funzione genera un vettore di numeri interi casuali Normali utilizzando il metodo Mersene-Twist.
%  

rng default

%tic

U1=rand(1,n);
U2=rand(1,n);

U=sqrt(-2*log(U1)).*cos(2*pi*U2);
clear U1 U2

%tempo_CONGL=toc ;
%disp(['Tempo calcolo: ' num2str(tempo_CONGL)])


end

