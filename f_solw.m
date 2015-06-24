function w = f_solw(Ls,wt,r)
% This function calculates equilibrium objects and calculates a new wage
% conditional on these equilibrium objects

global n s t L modsigma eta f phi mu tau total tole

%% Calculating many objects
% Calculating the prices (p_ijs)
[tradeshare,inc,secinc,secexp] = f_soleqobjects(Ls,wt,r);
%% Calculating new eq consistent wages

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

%Wage rate
w=tradeflows_agg./L;


