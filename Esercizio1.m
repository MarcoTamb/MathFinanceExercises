clear

tipo = 1; %1 Opzione con massimo, 2 Opzione con integrale.


r = 0.1;
sigma = 2;
S0 = 100;
K = 100;
T = 1;
t = 0;

N = 100000; %dimensione MC simulation
M = 500; %numbero di passi temporali. 


% Generazione moto browniano come in BSMC.m

    dW = randn(M,N);
    dt = (T-t)/M;
    S = S0*exp(cumsum((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*dW));
   
% 
%   S è una matrice NxM, in ciascuna delle N colonne, corrispondenti a N 
%   estrazioni indipendenti di un moto Browniano, abbiamo il valore 
%   dell'asset negli M tempi.
%
  
switch tipo
    case 1 % opzioni con massimo 
        fprintf('-------------------------------------------------\n');
        fprintf('OPZIONE CON MASSIMO              Valore opzione \n');
        fprintf('-------------------------------------------------\n');
         
        %Put_T = max(K-max(S),0);
        Call_T = max(max(S)-K,0);
       
       %
       % max di una matrice è un vettore (riga) con il massimo fatto per
       % colonna. 
       %
       % Il meno tra un vettore e uno scalare agisce elemento per
       % elemento.
       %
       % Il max tra un vettore e uno scalare agisce elemento per elemento.
       % 
       
    case 2 %opzioni con integrale
        
        fprintf('-------------------------------------------------\n');
        fprintf('OPZIONE CON INTEGRALE               Valore opzione \n');
        fprintf('-------------------------------------------------\n');
        
        %Put_T = max(K-trapz(dt,S),0);
        Call_T = max(trapz(dt,S)-K,0);
        
       %
       % trapz integra con il metodo del trapezio con spaziatura dt. trapz
       % di una matrice è un vettore (riga) con l'integrale fatto per
       % colonna.
       %
       % Il meno tra un vettore e uno scalare agisce elemento per
       % elemento.
       %
       % Il max tra un vettore e uno scalare agisce elemento per elemento.
       %        
end

%valori MC
CallMC = exp(-r*(T-t))*mean(Call_T);
%PutMC =  exp(-r*(T-t))*mean(Put_T);


%errore MC
MC_error = sqrt(exp(-r*(T-t)))*[std(Call_T)]./sqrt(N);



fprintf('approssimazione MC                        %5.4f \n',[CallMC]);
fprintf('stima errore MC (sigma/sqrt(n))             %5.4f \n',MC_error);
fprintf('-------------------------------------------------\n');