function [trade_flows,tradeflows_sec_store,sector,tausigmas,w_out,market_acc] = f_output(Ls,t)
% This function calculates equilibrium objects and puts them into a common
% matrix for output
global n s L modsigma eta f phi mu tau total tole distance rta
%% Generating empty tables,loop control variables and initiate loop
    trade_flows = zeros(eval('n*n*s*t'),5);
    sector = zeros(eval('n*s*t'),11);
    tausigmas = zeros(eval('n*n*s*t'),7);
    counterstore1=1;
    counterstore2=1;
    counterstore3=1;
    Lsstore=Ls;
for r=1:t,
    %% Calculating wages
    Ls=Lsstore(:,:,r);
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
%     wt =wt/wt(1);
%     wt(1)=1;
    w_out(:,r)=wt;

    %% Calculating many objects

    [tradeshare,inc,secinc,secexp,market_acc] = f_soleqobjects(Ls,wt,r);
    market_acc(:,:,r)=market_acc;
    % Different trade flows at different levels of aggregation
    secexpt=secexp';
    for i=1:s,
        tradeflows_sec(:,:,i)=eval('bsxfun(@times,tradeshare(:,:,i),secexpt(i,:))');
        tradeflows_sec_store(:,:,i,r)=tradeflows_sec(:,:,i);
    end

    tradeflows_sec_country=sum(tradeflows_sec,2);
    tradeflows_agg=sum(tradeflows_sec_country,3);
    tradeflows_sec_country_compact=squeeze(tradeflows_sec_country);

    
  
    %% Add measurement error to trade flows
    % Currently deactivitated
%     errorterm= normrnd(0,0.1,[n,n,s]);
%     
%     errorterm=ones(n,n,s);
%     for i=1:n,
%        errorterm(i,i,:)=0; 
%     end
% 
%     for i=1:s,
%         errorterm(:,:,i)=errorterm(:,:,i)./(1-modsigma(i));
%     end
%     
% 
%     errorterm=exp(errorterm);
%     for j=1:s,
%         errorterm(:,:,j)=errorterm(:,:,j)./tradeflows_sec(:,:,j);
%     end
%     for j=1:s,
%         tradeflows_sec(:,:,j)=tradeflows_sec(:,:,j).*(errorterm(:,:,j));
%     end
% 

    %% Outputting the data in sector table
    %Generate tradeflow table: [year, sector s, Origin country i, destination country j, trade flow T_ijs]
    counter=1;
    if counterstore1 > 1,
    	counter=counterstore1;
    end
    for i=1:n,
        for j=1:n,  
            for q=1:s,
                trade_flows(counter,1)=r;
                trade_flows(counter,2)=q;
                trade_flows(counter,3)=i;
                trade_flows(counter,4)=j;
                trade_flows(counter,5)=tradeflows_sec(i,j,q);
                counter=counter+1;
            end
        end
    end
    counterstore1=counter;
    
    %Generate trade cost table: [year, sector, origin country i,
    %destination country j, tausigma_ijs]
    
    minussigmat=transpose(1-modsigma);
    for i=1:s,
        tausigma(:,:,i)=eval('bsxfun(@power,tau(:,:,i,r),minussigmat(i))');
    end
    

    
    counter=1;
    if counterstore2 > 1,
    	counter=counterstore2;
    end
    for i=1:n,
        for j=1:n, 
            for q=1:s,
                tausigmas(counter,1)=r;
                tausigmas(counter,2)=q;
                tausigmas(counter,3)=i;
                tausigmas(counter,4)=j;
                tausigmas(counter,5)=tausigma(i,j,q);
                tausigmas(counter,6)=distance(i,j);
                tausigmas(counter,7)=rta(i,j,q,r);
                counter=counter+1;
            end
        end
    end
    counterstore2=counter;
    %Generate Sectoral data: (Country, sector, year, L, w, VA (=Phi*L))
    counter=1;
    if counterstore3 > 1,
    	counter=counterstore3;
    end
    va=phi(:,:,r).*Ls;
    for i=1:n, 
        for q=1:s,
            sector(counter,1)=i;
            sector(counter,2)=q;
            sector(counter,3)=r;
            sector(counter,4)=Ls(i,q);
            sector(counter,5)=wt(i);
            sector(counter,6)=va(i,q);
            sector(counter,7)=f(i,q);
            sector(counter,8)=mu(i,q);
            sector(counter,9)=phi(i,q,r);
            sector(counter,10)=L(i);
            sector(counter,11)=market_acc(i,q);
            counter=counter+1;
        end
    end
    counterstore3=counter;
    %% Ending loop over time period
end
