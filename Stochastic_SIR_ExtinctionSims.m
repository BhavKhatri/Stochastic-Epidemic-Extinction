function [t,S,I,R]=Stochastic_SIR_ExtinctionSims(I0,Rec0,Re,Gamma,N)

% Gillespie/Kinetic Monte Carlo SIR sims
% I0 initial number infected
% Rec0 initial number recovered
% N Total population size
% Gamma rate of recovery in (days)^-1
% Re the effective reproductive number (=β*S0/γ)


% %Simulations expect Re<=1
% if Re>1
%     disp('')
%     disp('Error: simulations require Re<=1')
%     
%     Re = input('Input value of Re<=1:');
% end


%Initial susceptible
S0 = N - I0-Rec0;

%"Rate" of growth reaction 
k1 = Gamma*Re/S0;
% k1 = Gamma*Re;
%Rate of decay reactiom
k2 = Gamma;

%Initial number of timesteps
M=1e6;

%Initialise arrays — these will grow if need be
t = zeros(1,M);
S = zeros(1,M);
I = zeros(1,M);
R = zeros(1,M);

r1 = rand(1,M);
r2 = rand(1,M);



t(1)=0;
tt=0;
%Initial infected
I(1)=I0; 

%Initial susceptible
S(1) = S0;
% S = S0*ones(size(S));

%Initial Recovered
R(1) = Rec0;

k=1;


%Run simulation until I=0
while I(k)~=0
    
   
     if k+1 > M
        t = [t,zeros(1,100)];
        I = [I,zeros(1,100)];
        S = [S,zeros(1,100)];
%         S = [S,S0*ones(1,100)];
        R = [R,zeros(1,100)];
        r1 = [r1,rand(1,100)];
        r2 = [r2,rand(1,100)];
        
        M=M+100;
     end
            
     

     
    w1 = k1*S(k)*I(k);
    w2 = k2*I(k);
    
   
   
    if r1(k)<= w1/(w1+w2)
        I(k+1) = I(k)+1;
        S(k+1) = S(k)-1;
        R(k+1) = R(k);
    else
        I(k+1) = I(k)-1;
        R(k+1) = R(k)+1;
        S(k+1) = S(k);
    end
    
    t(k+1) = t(k) + log(1/r2(k))/(w1+w2);

    
    if t(k+1)> tt
        
       disp(['t = ',num2str(t(k+1))])
       
    tt=tt+10;
        
    end
    
%     disp(['t = ',num2str(t(k+1))])

    k=k+1;

    
   
end





S = S(1:k);
I = I(1:k);
R = R(1:k);

t = t(1:k);




