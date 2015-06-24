function Ls = f_soleq(Ls,r)
%------------------------------------------------------------------%
% Monday, 28th April 2014
% Created by: Simon Fuchs
% This function calculates the general equilibrium of the augmented
% Krugmann Model given the exogenous parameters specified in the main file.

%% Loading variables and defining
global n s t L modsigma eta f phi mu tau total tole

%% Wage loop

wt=ones(n,1);
err = 1e+6;
iter = 1;
while err > tole
    if iter > 1
        wt = neww;
    end
    w = f_solw(Ls,wt,r);
    neww = 0.9*wt + 0.1*w;
    err = sqrt(sum((neww-wt).^2)); % L2 norm
    iter = iter + 1;        
    if iter > 10000
        display('Iteration max')
        break
    end
end

% Normalization 
% wt =wt/wt(1);
% wt(1)=1;
%     


%% Calculating many objects
[tradeshare,inc,secinc,secexp] = f_soleqobjects(Ls,wt,r);
%% Calculating new labor allocation

% Sectoral Goods/Labor market clearing
secexpt=secexp';


for i=1:s,
    tradeflows_sec(:,:,i)=eval('bsxfun(@times,tradeshare(:,:,i),secexpt(i,:))');
end

tradeflows_sec_country=sum(tradeflows_sec,2);

for i=1:s,
    market_wedge_sec_store(:,i)=secinc(:,i)-tradeflows_sec_country(:,i);
end

market_wedge_sec=(reshape(market_wedge_sec_store,1,[]))';

% National Goods/Labor market clearing
tradeflows_agg=sum(tradeflows_sec_country,3);
market_wedge_agg=inc-tradeflows_agg;

%Labor market clearing
labor_market_wedge=sum(Ls,2)-L;

tradeflows_sec_country_compact=squeeze(tradeflows_sec_country);

for i=1:s,
    Ls(:,i)=eval('bsxfun(@ldivide,wt,tradeflows_sec_country_compact(:,i))');
end


