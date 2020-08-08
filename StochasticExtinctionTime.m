<<<<<<< HEAD
function [tau,sig_tau,medtau,ltau,utau,CItau] = StochasticExtinctionTime(Re,Gamma,I0,S0,CI, plotfig)
=======
function [tau,sig_tau,medtau,ltau,utau,CItau] = StochasticExtinctionTime(Re,Gamma,I0,CI, plotfig)
>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043

%function that returns the mean (tau), standard devation (sig_tau), median
%(medtau) and confidence intervals (ltau, utau) as specified by 0<CI<1, for
%each value of Re given in vector/array Re

%To plot the distributions set plotfig=1 — if the number of elements of Re
%> 5 then plot is not peformed as there would be too many figures on one
%graph. Note that if 1-Re is very small extinction time will be very large
%compared to say Re=0.5 — so setting the Xscale to log might help
% e.g. set(gca,'XScale','log')

%Check
if Re>=1
    fprintf('\n')
    disp('**#!Error to plot distribution of times to extinction Re<1 — Re>1 gives growth of infections!')
    fprintf('\n')
    return
elseif Re<0
    fprintf('\n')
    disp('Error Re cannot be negative!')
    fprintf('\n')
    return
end

rho = Gamma*(1-Re);

Idagger = 1./(1-Re);
tdagger = 1./rho.*log(I0./Idagger);

%Mean and standard deviation
tau = double(eulergamma)./rho + tdagger;

sig_tau = pi/sqrt(6)./rho;

%Median
medtau = tdagger - 1./rho*log(log(2));


%Confidence Interval
alpha = 1/2*(1-CI);

ltau = tdagger - 1./rho*log(-log(alpha));
utau = tdagger - 1./rho*log(-log(1-alpha));


if numel(Re)==1
    CItau = [ltau,medtau,utau];
else
    CItau = [];
end

<<<<<<< HEAD
l=1;

h = [];
=======

>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043
if numel(Re)<=5 && plotfig==1
    
    figure; 
    
    for n=1:numel(Re)
        
        if numel(Re)==1
        
<<<<<<< HEAD
            H=plotdist(rho(n),Re(n),tdagger(n),I0,Idagger(n),tau(n),medtau(n),ltau(n),utau(n),alpha,n,1); 
            
            
            if ~isempty(S0) && ( 1-Re < 2*sqrt(I0/S0) )

                fprintf('\n')
                disp(['Warning: for Re =',num2str(Re(n))])
                disp('Re may be too close to 1 for original approximation of unchanging susceptible population to be valid')
            end
            
            
            if numel(H)==1
                legendstr{l} = ['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho,2)),'\ \mathrm{\%reduction/day}$'];
                h = [h,H];
                               
                l=l+1;
            else
                legendstr{l}=['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho,2)),'\ \mathrm{\%reduction/day}$ (Approx)'];
                legendstr{l+1}=['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho,2)),'\ \mathrm{\%reduction/day}$ (Exact)'];
                
                h = [h,H];
                
                l=l+2;
            end
        
        else
            
                       
            H=plotdist(rho(n),Re(n),tdagger(n),I0,Idagger(n),tau(n),medtau(n),ltau(n),utau(n),alpha,n,0);
            
            if ~isempty(S0) & (1-Re(n) < 2*sqrt(I0/S0) )

                fprintf('\n')
                disp(['Warning: for Re = ',num2str(Re(n))])
                disp('Re may be too close to 1 for original approximation of unchanging susceptible population to be valid')
                disp('    Estimation of extinction times likely to be an overestimate')
            end
            
            if numel(H) == 1
            
                legendstr{l} = ['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re(n)),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho(n),2)),'\mathrm{\%\ reduction/day}$'];
                
                h = [h,H];
                               
                l=l+1;
                
            else
                legendstr{l}=['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re(n)),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho(n),2)),'\ \mathrm{\%reduction/day}$ (Approx)'];
                legendstr{l+1}=['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re(n)),...
                    ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
                    '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho(n),2)),'\ \mathrm{\%reduction/day}$ (Exact)'];
                
                h = [h,H];
                               
                l=l+2;
                
            end
                
=======
            H(n)=plotdist(rho(n),tdagger(n),tau(n),medtau(n),ltau(n),utau(n),alpha,n,1);            
            
            legendstr{n}=['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re),...
            ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
            '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho,2)),'\ \mathrm{\%reduction/day}$'];
        
        else
            
            H(n)=plotdist(rho(n),tdagger(n),tau(n),medtau(n),ltau(n),utau(n),alpha,n,0);
            
            legendstr{n} = ['$\mathrm{Reproductive\ number}\ R_e = ',num2str(Re(n)),...
            ';\ \ \mathrm{infection\ duration}\ 1/\gamma = ',num2str(1/Gamma),...
            '\ \mathrm{days};\quad \rho_e = ',num2str(round(100*rho(n),2)),'\mathrm{\%\ reduction/day}$'];
>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043
            
        end
        
        
    end
    
    
        
   
    xlabel 'Extinction time (days)'
    ylabel 'Probability density'
    
<<<<<<< HEAD
    legend(h,legendstr)
    
    title(['Stochastic SIR model (constant $R_e$ assumption) --- $1/\gamma = ',...
        num2str(round(1/Gamma,3)),'\ \mathrm{days}; I_0 = ',num2strpow(I0),'; S_0 = ',num2strpow(S0),'$'])

    
    
elseif numel(Re)>5 && plotfig==1
        
    disp('Too many values of Re to plot pdf on one plot')
=======
    legend(H,legendstr)
    
    elseif numel(Re)>5 && plotfig==1
            disp('Too many values of Re to plot pdf on one plot')
>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043
            
    
end


         
    
    
    
    

    
    
    
end



    
<<<<<<< HEAD
function h = plotdist(rho,Re,tdagger,I0,Idagger,tau,medtau,ltau,utau,alpha,n,addtext)
=======
function h = plotdist(rho,tdagger,tau,medtau,ltau,utau,alpha,n,addtext)
>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043



    load DefaultColourOrder.mat
    
    
    t=0:0.1:utau*1.5;
    p = rho*exp(-rho*(t-tdagger)).*exp(-exp(-rho*(t-tdagger)));
    
    h=plot(t,p,'-','Color',defcolours(n,:),'LineWidth',2);hold on
<<<<<<< HEAD
    
    if I0/Idagger>25
        
        col= defcolours(n,:);
        
        t=0:0.1:utau*1.5;
        p = rho*exp(-rho*(t-tdagger)).*exp(-exp(-rho*(t-tdagger)));
    
        h=plot(t,p,'-','Color',col,'LineWidth',2);hold on
        
        
        pp1 = interp1(t,p,tau);
        line([tau, tau],[0 pp1],'Color',col,'LineStyle','--','LineWidth',2)


        pp2 = interp1(t,p,medtau);
        line([medtau, medtau],[0 pp2],'Color',col,'LineStyle','--','LineWidth',1)


        pp3 = interp1(t,p,ltau);
        line([ltau, ltau],[0 pp3],'Color',col,'LineStyle','--','LineWidth',0.5)


        pp4 = interp1(t,p,utau);
        line([utau, utau],[0 pp4],'Color',col,'LineStyle','--','LineWidth',0.5)

        if addtext==1

            text(tau,1.025*pp1,'Mean','HorizontalAlignment','left')
            text(medtau,1.025*pp2,'Median','HorizontalAlignment','left')
            text(0.99*ltau,1.025*pp3,[num2str(100*alpha),'th \%tile'],'HorizontalAlignment','right')
            text(utau,1.2*pp4,[num2str(100*(1-alpha)),'th \%tile'],'HorizontalAlignment','left')

        end
    
    
        
    else
        
        fprintf('\n')
        disp(['Warning: for Re = ',num2str(Re)])
        disp('Gumbel approximation could be poor as I0 ~ Idagger (Idagger =1/(1-Re))')
        disp('Exact distribution plotted in dashed lines')
        disp('    — use judgement to assess whether mean, standard deviation, median and confidence intervals are sensible')
        disp('Mean of both are indicated')
        
        t=0:0.1:utau*1.5;
        p = rho*exp(-rho*(t-tdagger)).*exp(-exp(-rho*(t-tdagger)));
    
        h(1)=plot(t,p,'-','Color',defcolours(n,:),'LineWidth',2);hold on
        
        pp1 = interp1(t,p,tau);
        line([tau, tau],[0 pp1],'Color',defcolours(n,:),'LineStyle','--','LineWidth',1)
        
        
        
        tq=0:0.1:utau*10;
        q = (1-Re)*rho*I0*exp(rho*tq).*(1+(Re-1)./(exp(rho*tq)-Re)).^(I0-1)./(exp(rho*tq)-1).^2;
        
        h(2) = plot(tq,q,'--','Color',defcolours(n,:),'LineWidth',2);
        
        ind = ~isnan(q);
        tauq = sum(tq(ind).*q(ind))*0.1;
        pp2 = interp1(tq,q,tauq);
        line([tauq, tauq],[0 pp2],'Color',defcolours(n,:),'LineStyle',':','LineWidth',1)
        
        xlim([0, utau*1.5])
        
        
        
        
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

=======
        
    
    hfill=fill(t,p,defcolours(n,:));
    set(hfill,'FaceAlpha',0.4)
    
%     u = axis;
    
    pp1 = interp1(t,p,tau);
    line([tau, tau],[0 pp1],'Color','k','LineStyle','--','LineWidth',2)
    
    
    pp2 = interp1(t,p,medtau);
    line([medtau, medtau],[0 pp2],'Color','k','LineStyle','--','LineWidth',1)
    
    
    pp3 = interp1(t,p,ltau);
    line([ltau, ltau],[0 pp3],'Color','k','LineStyle','--','LineWidth',0.5)
    
    
    pp4 = interp1(t,p,utau);
    line([utau, utau],[0 pp4],'Color','k','LineStyle','--','LineWidth',0.5)
    
    if addtext==1
      
        text(tau,1.025*pp1,'Mean','HorizontalAlignment','left')
        text(medtau,1.025*pp2,'Median','HorizontalAlignment','left')
        text(0.99*ltau,1.025*pp3,[num2str(100*alpha),'th \%tile'],'HorizontalAlignment','right')
        text(utau,1.2*pp4,[num2str(100*(1-alpha)),'th \%tile'],'HorizontalAlignment','left')
    
        
        
    
    
    
    end
    
>>>>>>> 6589ae5ead2213f70678e44e6c32869f5bd2d043

end

