%Valutazione Monte Carlo di opzioni
clear

tipo = 0; %1 Esotica doppio massimo, 2 Esotica integrale, 0 Europea. 


r = 0.1;
sigma = 2;
T = 1;
t = 0;
S0 = 100;
K = 100;


N = 1000000; %dimensione MC simulation
M = 252; %numbero di passi temporali. 


% Generazione moto browniano come in BSMC.m
    dW = randn(M,N);
    dt = (T-t)/M;
    S = S0*exp(cumsum((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*dW));
  
  % S è una matrice NxM, in ciascuna delle N colonne abbiamo il valore 
  % dell'asset negli M tempi.
  
switch tipo
    case 1 % opzioni esotiche "compito per casa" con massimo 
        
       fprintf('OPZIONE CON MASSIMO \n');
       Put_T = max(K-max(S),0);
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
       
    case 2 %opzioni esotiche "compito per casa" con integrale
        
        fprintf('OPZIONE CON INTEGRALE \n');
        
        Put_T = max(K-trapz(dt,S),0);
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
        
    case 0 %opzioni europee
        fprintf('OPZIONE EUROPEA\n');
        
        Put_T = max(K-S(M,:),0);
        Call_T = max(S(M,:)-K,0);
        
        %
        % Il meno tra un vettore e uno scalare agisce elemento per
        % elemento.
        %
        % Il max tra un vettore e uno scalare agisce elemento per elemento.
        %
        
end
%valori MC
CallMC = exp(-r*(T-t))*mean(Call_T);
PutMC =  exp(-r*(T-t))*mean(Put_T);

%errore MC
MC_error = sqrt(exp(-r*(T-t)))*[std(Call_T) std(Put_T)]./sqrt(N);

%[Call,Put] = blsprice(S0,K,r,T-t,sigma); %Valore BS esatto

fprintf('-----------------------------------------------\n');
fprintf('                                 Call   Put\n');
fprintf('-----------------------------------------------\n');
%fprintf('valore esatto opz. europee     %5.4f %5.4f\n',[Call Put]);
fprintf('approssimazione MC             %5.4f %5.4f\n',[CallMC PutMC]);
fprintf('errore MC                      %5.4f %5.4f\n',MC_error);
fprintf('-----------------------------------------------\n');