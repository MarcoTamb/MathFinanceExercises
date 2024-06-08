function [U] = PRNGNormalGraph(n)
%PRNGUNIFORMCONG Questa funzione genera un vettore di numeri interi casuali Normali utilizzando il metodo Mersene-Twist.
%  

rng default

tic

U1=rand(1,n);
U2=rand(1,n);

U=sqrt(-2*log(U1)).*cos(2*pi*U2);
clear U1 U2

tempo_CONGL=toc ;
disp(['Tempo calcolo: ' num2str(tempo_CONGL)])

% Grafici
figure
plot(U(1:2:end-1),U(2:2:end),'.')
title('Numeri pseudocasuali - Normale')
xlabel('x_{2i}')
ylabel('x_{2i+1}')

figure
histogram(U,50)
title('Normale - Istogramma')

figure
cdfplot(U)
grid
xlabel('u.i.d.')
title('Normale - Funzione di ripartizione')

end

