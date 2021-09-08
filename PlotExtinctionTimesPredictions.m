function PlotExtinctionTimesPredictions(Gamma,I0,S0,Re,N,Rec0,CI,col)


if isempty(Re)
    Re = [0:0.005:1-0.00005];
    
end

if isempty(col)
    col = [0.466        0.674        0.188];
end

if isempty(CI)
    CI = 0.95;
end

while  isempty(S0)
    disp('Value of initial susceptible population size required to assess validity of predictions')
    S0=input('Please input value:')
end


%Calculate extinction times
[tau,sig_tau,medtau,ltau,utau,~] = StochasticExtinctionTime(Re,Gamma,I0,S0,CI,0,0);

figure;

%Plot mean and confidence intervals
[h0,h1]=plot_errorbound(Re,tau,tau-ltau,utau-tau,col,0);hold on

ymax = 1000;
ylim([0 ymax])


%Plot median
hmed = plot(Re,medtau,'--','Color',col);


%Plot deterministic prediction
rho = Gamma*(1-Re);
tau_det = 1./rho*log(I0);
hdet = plot(Re,tau_det,'.','Color',col);



%Re close to 1 correction
x=0.95;
% Restar = S0*log(x)/(Rec0+S0*x -N)




% re = Restar:0.005:1;
re = Re;
r0 = N*re/S0;

ReInf = -lambertw(-re.*exp(-r0.*(1-Rec0/N)));
rhoeInf = Gamma*(1-ReInf);
% Restar = ReInf(end)
Restar = -lambertw(-exp(-N/S0.*(1-Rec0/N)))

% if max(ReInf./Re)<x
%     x= max(ReInf./Re)*1.05
%     Restar = log(1/x)/(N/S0*(1-Rec0/N)-x)
    
    
% Restar = 1-sqrt(I0/S0);

[tau1,sig_tau1,medtau1,ltau1,utau1,~] = StochasticExtinctionTime(ReInf,Gamma,I0,S0,CI,0,0);
% tau1

   
% indtau = re>=Restar;
% indtau = rhoeInf<1./tau

% tau1(indtau) = tau1(indtau)+1./rhoeInf(indtau);
% ltau1(indtau) =ltau1(indtau)+1./rhoeInf(indtau);
% utau1(indtau)  = utau1(indtau)+1./rhoeInf(indtau);


sw=switchfunc(Re,Restar,0.05);

Dtau = sw./rhoeInf;

tau1 = tau1+Dtau;
ltau1 =ltau1+Dtau;
utau1  = utau1+Dtau;

% tau1 = tau1+1./rhoeInf;
% ltau1 =ltau1+1./rhoeInf;
% utau1  = utau1+1./rhoeInf;


hRe1 = plot(re,tau1,'-.','Color',col);

%Critical Re corresponding to when I0 = 25*Idagger, where Idagger=1/(1-Re)
%is the stochastic threshold
RRe = 1-25/I0;

r = [RRe, RRe,1,1];

hnotexact = fill(r,[0 ymax ymax 0],[0.8 0.8 0.8]);
set(hnotexact,'FaceAlpha',0.4)

    
%Plot exact mean from branching process theory with no approximations (mean
%calculated numerically in function below from exact distribution
for n=1:numel(Re)
    
    tau_exact(n) = MeanExactExtinctionTimeDistribution(Re(n),Gamma,I0);
    
end

hexact = plot(Re,tau_exact,'k','LineWidth',1);


%Plot vertical line to indicate value of Re that constant S(t)
%approximation breaks down
% Restar = 1-sqrt(I0/S0);
hsmallRe = line([Restar, Restar],[0,ymax],'Color','k','LineStyle','--','LineWidth',1);





%Legend
h = [h1,hmed,h0, hexact,hdet,hRe1,hnotexact,hsmallRe];

legendstr{1} = 'Mean (Gumbel)';
legendstr{2} = 'Median (Gumbel)';
legendstr{3} = '95\% CI (Gumbel)';
legendstr{4} = 'Mean (Exact)';
legendstr{5} = 'Deterministic';
legendstr{6} = 'Mean Time Correction for $R_e>\sqrt{I_0/S_0}$';
legendstr{7} = 'Region: $I_0\sim I^\dagger\ \Rightarrow$ Gumbel distribution not exact';
legendstr{8} = 'Small $1-R_e$ threshold: When $1-R_e\sim < \sqrt{I_0/S_0}$ constant $R_e$ approximation poor';
    
legend(h,legendstr,'Location','northwest');


title(['Stochastic SIR model (no herd immunity) --- $1/\gamma = ',...
        num2str(round(1/Gamma,3)),'\ \mathrm{days}; I_0 = ',num2strpow(I0),'; S_0 = ',num2strpow(S0),'$'])





hold off
grid on

xlabel('Effective Reproductive Number $R_e$')
ylabel('Extinction Time (days)')


ind1 =tau<=tau1;
ind2 = tau>tau1;

Tau = [tau(ind1),tau1(ind2)];
sig_Tau  = [sig_tau(ind1),sig_tau1(ind2)];
lTau = [ltau(ind1),ltau1(ind2)];
uTau = [utau(ind1),utau1(ind2)];

% sw=switchfunc(Re,Restar,0.05)
% 
% Tau = tau.*(1-sw) + tau1.*sw;
% lTau = ltau.*(1-sw) + ltau1.*sw;
% uTau = utau.*(1-sw) + utau1.*sw;

figure;

%Plot mean and confidence intervals
[H0,H1]=plot_errorbound(Re,Tau,Tau-lTau,uTau-Tau,col,0);%hold on

xlabel('Effective Reproductive Number $R_e$')
ylabel('Extinction Time (days)')
title(['Stochastic SIR model (no herd immunity) --- $1/\gamma = ',...
        num2str(round(1/Gamma,3)),'\ \mathrm{days}; I_0 = ',num2strpow(I0),'; S_0 = ',num2strpow(S0),'$'])

% indtau = find(abs(tau-tau1)<0.01);
% RRestar = mean

figure;plot(Re,ReInf)
end


function f=switchfunc(x,a,b)

f=0.5 +0.5*tanh((x-a)/b);

end



function [h0,h1]=plot_errorbound(x,y,ep,em,col,bound_only)


    
    ep = y+ep;
    em = y-em;

    ep=squeeze(ep);
    em=squeeze(em);
%     size(em)

    indNotNan = ~isnan(em);
    
    em = em(indNotNan);
    ep = ep(indNotNan);
    xx = x(indNotNan);
    
    %Create vector to fill graph

    if isrow(em)
        e=[em,fliplr(ep)];
    else
        e=[em;flipud(ep)];
    end

    if isrow(xx)
        xx=[xx,fliplr(xx)];
    else
        xx=[xx;flipud(xx)];
    end



    h0=fill(xx,e,col);hold on
    set(h0,'FaceAlpha',0.2,'LineStyle','none')
    
    if bound_only==0
        h1=plot(x,y,'LineWidth',2,'LineStyle','-');%hold on
        set(h1,'Color',col); %hold off
    end
    


end


function S = num2strpow(s)

pow = floor(log10(s));

if pow <5
    S = num2str(s);
else
    C = round(10^(log10(s)-pow),3);

    S = [num2str(C),'\times',' 10^',num2str(pow)];
    
end


end





