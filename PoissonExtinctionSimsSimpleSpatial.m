function [t,I]=PoissonExtinctionSimsSimpleSpatial(I0,Re,gamma,m,n)

%Gillespie/Kinetic Monte Carlo sims of n subpopulations with global migration
%probability m

% I0 initial number infected
% Gamma rate of recovery in (days)^-1
% Re the effective reproductive number

%Simulations expect Re<1
if Re>1
    disp('')
    disp('Error: simulations require Re<=1')
    
    Re = input('Input value of Re<=1:');
end




%Rate of growth reaction
k1 = gamma*Re;
%Rate of decay reactiom
k2= gamma;

% %Arrays/vectors of rates for each sub-population
% kplus = zeros(1,n);
% kminus = zeros(1,n);
% k = zeros(1,2*n);

%Max #timesteps
Nt=1e6;



t = zeros(1,Nt);
I = zeros(n,Nt);

%Migration probability for each sub-population
mout = m*ones(n,1);

%Initialise migration matrix
M = zeros(n);

%populate migration matrix — here assumption is global migration
%(probability of migration out of sub-population is m, probability to
%migrate to any other sub-population is m/(n-1)
for i=1:n
    for j=1:n
        if i~=j
            M(i,j) = mout(j)/(n-1);
        end
    end
end

M



t(1)=0;

%Assume initial infected are uniformly spread over sub-populations
I(:,1)=round(I0/n)*ones(n,1); 


maxsteps = Nt;

tk=1;

while sum(I(:,tk))~=0 %Condition for extinction is that no infected in any sub-population
    
   
     if tk+1 > maxsteps
        t = [t,zeros(1,Nt)];
        I = [I,zeros(n,Nt)];
        maxsteps = maxsteps+Nt;
     end
     
     II = I(:,tk);     

%      size(M)
%      size(II)
%      size(mout)
     
     kplus = k1*II;
     kminus = k2*II;
     
     kmplus = M*II;
     kmminus = mout.*II;
     
     k = [kplus;kminus;kmplus;kmminus];
     K=cumsum(k);
     %Find index which is closest to rand*totalrate
     u=rand;
     Kind=find(K>u*max(K),1); %first value greater excludes the zero rate wt states
     
     if Kind <= n
         %I-> I+1 from infection
%          Kind = Kind;
         I(:,tk+1) = I(:,tk);
         I(Kind,tk+1) = I(Kind,tk)+1;
     elseif Kind > n & Kind <= 2*n
         %I-> I-1 from recovery
         Kind = Kind-n;
         I(:,tk+1) = I(:,tk);
         I(Kind,tk+1) = I(Kind,tk)-1;
     elseif Kind > 2*n & Kind <= 3*n
         %I-> I+1 from emigration
         Kind = Kind-2*n;
         I(:,tk+1) = I(:,tk);
         I(Kind,tk+1) = I(Kind,tk)+1;
     else
         %I-> I-1 from immigration
         Kind = Kind-3*n;
         I(:,tk+1) = I(:,tk);
         I(Kind,tk+1) = I(Kind,tk)-1;
     end
  
    
    t(tk+1) = t(tk) + log(1/rand)/max(K);

    
    tk = tk+1;
  
   
    
end

I = I(:,1:tk);
t = t(1:tk);




