clear 

sigma=2;
T=1;
R=0.1; %r=RT/N, S^0_n=(1+r)^n.

S_0=100; %valore iniziale
K=100; %prezzo di esercizio

% Impostazione parametri 
N = 1000; %numbero di passi temporali. 
r=R*T/N; %S^0_n=(1+r)^n, r deve essere in (d-1, u-1)
u=exp(sigma*sqrt(T/N))*(1+r); %incremento  
d=exp(-sigma*sqrt(T/N))*(1+r); %decremento

% Per eseguire come simulazione montecarlo con u, d e r generici, 
% commentare r, u e d  qui sopra e inserirli esplicitamente qui sotto, 
% rimuovendo il commento.
%
% r=
% u=
% d=
%
% if d>=u % necessario u>d
%     error('è necessario che valga u>d')
% end
% if r>(u-1)
%     error('r>u-1. Presenza di arbitraggio.')
% elseif r<(d-1)
%     error('r>d-1. Presenza di arbitraggio.')
% end
%


MC = 1000000; %dimensione MC simulation

% STIMA CON CDF NORMALE

sst=sigma.*sqrt(T);
t_1=(log(S_0/K)+R*T+(0.5*(sigma.^2).*T))./(sst);
t_2=t_1 - sst;
phi1=normcdf(t_1, 0, 1);
phi2=normcdf(t_2, 0, 1);
CN=S_0*phi1-(K)*exp(-R*T)*phi2;

% MONTECARLO Cox-Ross-Rubinstein

%generazione della matrice contenente gli incrementi

p=1-((u-1-r)/(u-d));        %probability of u

Tt=rand([N,MC]); %N righe da T elementi

Tt=u*(Tt<p)+d*(Tt>p); %Tt(i, j) è l'incremento all'i-esimo tempo nella j-esima simulazione
S=S_0*prod(Tt);

clear Tt; % clear per Tt perché usa troppa RAM.

% Valore al tempo T di h
hT=max((S-K), 0);


% Valore scontato medio
h=mean(hT)*((1+r)^-N);

%stima errore
MC_error = sqrt((1+r)^-N)*std(hT)./sqrt(MC);

fprintf('----------------------------------------------------------\n');
fprintf('                                               Valore Call \n');
fprintf('----------------------------------------------------------\n');
fprintf('approssimazione Call con formula C-R-R:           %5.4f \n',[CN]);
fprintf('approssimazione MC C-R-R con N=%d passi          %5.4f \n',[N h]);
fprintf('stima errore MC                                    %5.4f\n',MC_error);
fprintf('Numero iterazioni MC                              %d \n',[MC]);
fprintf('----------------------------------------------------------\n');
clear d hT K MC N p phi1 phi2 r R S S_0 sigma sst T t_1 t_2 u