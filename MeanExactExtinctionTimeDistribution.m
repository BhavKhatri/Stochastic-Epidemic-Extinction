function tau = MeanExactExtinctionTimeDistribution(Re,Gamma,I0)

        [~,~,~,~,utau,~] = StochasticExtinctionTime(Re,Gamma,I0,[],0.95,0);
        
        rho = Gamma*(1-Re);
        
        t=0:0.1:utau*10;
        p = (1-Re)*rho*I0*exp(rho*t).*(1+(Re-1)./(exp(rho*t)-Re)).^(I0-1)./(exp(rho*t)-1).^2;
        
%         h(2) = plot(t,q,'--','Color',defcolours(n,:),'LineWidth',2);
        
        ind = ~isnan(p);
        tau = sum(t(ind).*p(ind))*0.1;
        
        
end


