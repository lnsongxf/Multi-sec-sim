%------------------------------------------------------------------%
% Monday, 05th November 2014
% Created by: Simon Fuchs
% 
% This programme simulates the model to generate a test dataset for 
% evaluating the estimation algorithm.
% 
%% Prelims
clc
clear 
cd('/Users/Lucks/Desktop/RA Simulation 07-11')
%% Exogenous parameters
global n s t L modsigma eta f phi mu tau total tole distance rta
%Convergence for internal loops 
tole=1e-13;

%Setting general parameters
n=25; %Number of countries
s=2; %Number of sectors
t=2; %Number of years

%Setting country specific parameters
L=normrnd(50,50,[n,1]);
L=abs(L)

%Setting sectoral global parameters
modsigma=rand(s,1)*4+1;
modsigma=[5.5 2.5]'
eta=rand(s,1)+0.2;
eta=[0.7 0.3]'
while max(-modsigma)>-2,
    modsigma=rand(s,1)*4+1;
end
modsigma
while max(eta)>1,
    eta=rand(s,1)+0.2;
end
eta
%Setting sectoral national parameters (for all i, for all s) 
% f=5*ones(n,s)
% f=normrnd(5,4,[n,s]);
% f=abs(f)
phi=normrnd(3,1,[n,s]);
phi=abs(phi)+1;
f=1.*(1./phi);


errorterm2= normrnd(0,0.2,[n,s,t]);
errorterm2(:,:,1)= ones([n,s]);
shift=(1+rand(t,1)*2);
shift(1)=1;
%errorterm2=ones([n,s,t]);
% for i=2:t,
% errorterm2(:,:,i)=errorterm2(:,:,i)*shift(i);
% end
phi2=ones([n,s,t])
for i=1:t,
    phi2(:,:,i)=phi+errorterm2(:,:,i);
end

phi=phi2

clear phi2;
mu=(ones(s,n)*1/s)';

%Generating trade costs from distance and regional trade agreements
%characteristics + disturbance term

%Case 1: Linear trade costs
% distance=normrnd(1.7,0.1,[n,n]);
% for i=1:n,
%     distance(i,i)=1;
% end
% distance=abs(distance);
% rta=normrnd(1,0.1,[n,n,s,t]);
% for i=1:n,
%     rta(i,i,:,:)=1;
% end
% tau=ones(n,n,s,t)*0.9;
% while min(min(min(min(tau))))<1,
%     beta1=1
%     for j=1:s,
%         for i=1:t,
%             tau(:,:,j,i)=beta1*distance+rta(:,:,j,i);
%         end
%     end
%     min(min(min(min(tau))))
% end

%Case 2: Exponential trade costs
% 
distance=normrnd(1,0.1,[n,n]);
% for i=1:n,
%     distance(i,i)=0;
% end
distance=abs(distance);
rta=normrnd(0,0.2,[n,n,s,t]);
rta=zeros([n,n,s,t]);
% for i=1:n,
%     rta(i,i,:,:)=0;
% end
for i=1:n,
    rta(:,:,:,1)=0;
end
tau=ones(n,n,s,t)*0.9;

beta2=1./(1-modsigma);

for i=1:s,
    rta2(:,:,i,:)=(1+rta(:,:,i,:)).^beta2(i);
end
shift2=(1+0.1*rand(t,1));
shift2(1)=1;
%shift2=ones([1,t]);

while min(min(min(min(tau))))<1,
    beta1=0.55.*shift2;
    beta3=1.7;
    for j=1:s,
        for i=1:t,
            tau(:,:,j,i)=exp(beta1(i)*distance);
        end
    end
%     for j=1:s,
%         for i=1:t,
%             tau(:,:,j,i)=beta3*distance;
%         end
%     end
%    for j=1:n,
%             tau(j,j,:,:)=1;
%     end
%     for j=1:s,
%         for i=1:t,
%             rta(:,:,j,i)=rta(:,:,j,i)./tau(:,:,j,i);
%         end
%     end
    for j=1:s,
        for i=1:t,
%             tau(:,:,j,i)=tau(:,:,j,i)+rta(:,:,j,i);
            tau(:,:,j,i)=tau(:,:,j,i).*(rta2(:,:,j,i));
        end
    end
    min(min(min(min(tau))))
end



%% Solving for equilibrium wage/labor allocation
% This section iterates on equilibrium conditions. The outer loop exploits
% sectoral clearing to solve for the labor allocation across sectors; the
% inner loop uses aggregate clearing condition to solve for the eq
% consistent wage level. 
for r=1:t,
    % Set initial labor allocation to aggregate labor evenly distributed
    total=n*s;
    x = ones(total,1)';  
    x(eval('1'):eval('s'))=L(1)/s;
    for i=2:n,
       x(eval('1+(i-1)*s'):eval('i*s'))=L(i)/s;
    end
    Ls_old=[x(eval('1'):eval('s'))];
    for i=2:n,
       Ls_old=[Ls_old; x(eval('1+(i-1)*s'):eval('i*s'))];
    end
    clear total x;
    
    %Set iteration parameters
    err = 1e+6;
    iter = 1;
    
    %Initiate loop
    while err > tole
        if iter > 1
            Ls_old = Ls_new;
        end
        Ls = f_soleq(Ls_old,r);
        Ls_new = 0.6*Ls_old + 0.4*Ls;
        err = max(max(sqrt((Ls_new-Ls_old).^2))) % L2 norm
        iter = iter + 1;        
        if iter > 10000
            display('Iteration max')
            break
        end
    end
    Ls_out(:,:,r)=Ls;
end
Ls=Ls_out;
%% Outputting the data
[trade_flows,tradeflows_sec,sector,tausigmas,w,market_acc] = f_output(Ls,t);
[tausigmas_err]=f_output_error(Ls,t);
cd Output/
headers = {'year','sitc','origin','destination','value'};
csvwrite_with_headers('tradeflows.csv',trade_flows,headers);

headers = {'country','sitc','year','employment','wage','value added','fix cost','exp share','productivity','country employment','market access'};
csvwrite_with_headers('sector.csv',sector,headers);

% Generate Simulation data
headers = {'country no','sector numbers','year'};
simbasic=[n,s,t];
csvwrite_with_headers('sim.csv',simbasic,headers);
L_out=[1:n;L']'
headers = {'year','sitc','origin','destination','tausigma','distance','rta'};
csvwrite_with_headers('tausigma.csv',tausigmas,headers);
headers = {'year','sitc','origin','destination','tausigma'};
csvwrite_with_headers('tausigma-err.csv',tausigmas_err,headers);
csvwrite('fix.csv',f);
csvwrite('phi.csv',phi);
csvwrite('L.csv',L_out);
csvwrite('Ls.csv',Ls);
csvwrite('taus.csv',tau);
csvwrite('mu.csv',mu);
csvwrite('sigma.csv',modsigma);
modgamma=(1-modsigma).*beta1;
csvwrite('gamma.csv',modgamma);
csvwrite('eta.csv',eta);

%% Some tests concerning the wage equation approach.

tradeflows_sec_country=sum(tradeflows_sec,2);
tradeflows_sec_country_compact=squeeze(tradeflows_sec_country);
sector(1,5)*sector(1,4)-tradeflows_sec_country_compact(1,1,1)

LHS=(sector(1,5)/sector(3,5))^modsigma(1);

RHS2=(sector(3,7)/sector(1,7))*(sector(3,4)/sector(1,4))^eta(1)*(sector(1,9)/sector(3,9))^(1-modsigma(1))*(sector(1,11)/sector(3,11));

RHS=(sector(3,4)/sector(1,4))^eta(1)*(sector(1,9)/sector(3,9))^modsigma(1)*(sector(1,11)/sector(3,11));

LHS-RHS
LHS-RHS2


for i=1:s,
    wpower(:,i)=w(:,1).^modsigma(i);
    rel_lab(i)=(Ls(5,i,1)./Ls(3,i,1)).^eta(i);
    rel_phi(i)=(phi(5,i,1)./phi(3,i,1)).^(1-modsigma(i));
    rel_phi_adj(i)=(phi(3,i,1)./phi(5,i,1)).^(modsigma(i));
end

LHS=wpower(3,:)./wpower(5,:)
RHS=rel_lab.*(f(5,:)./f(3,:)).*rel_phi.*(market_acc(3,:,1)./market_acc(5,:,1))
RHS=rel_lab.*rel_phi_adj.*(market_acc(3,:,1)./market_acc(5,:,1))

