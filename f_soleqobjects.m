function [tradeshare,inc,secinc,secexp,market_acc] = f_soleqobjects(Ls,wt,r)
% This function calculates equilibrium objects and calculates a new wage
% conditional on these equilibrium objects

global n s t L modsigma eta f phi mu tau total tole

%% Calculating many objects
% Calculating the prices (p_ijs)

etat=eta';
unitcost=eval('bsxfun(@ldivide,phi(:,:,r),wt)');
markup=(modsigma./(modsigma-1));
markupt=markup';
for i=1:s,
    nonmuprice(:,:,i)=eval('bsxfun(@times,tau(:,:,i,r),unitcost(:,i))');
end
for i=1:s,
    pricing(:,:,i)=eval('bsxfun(@times,markupt(i),nonmuprice(:,:,i))');
end
minussigmat=transpose(1-modsigma);
for i=1:s,
    pricingsigma(:,:,i)=eval('bsxfun(@power,pricing(:,:,i),minussigmat(i))');
end
% Calculating size of sector/number of firms (M_is)
% -------------------------------- checked -------------------------------%

secsize=eval('bsxfun(@power,Ls,1-etat)')./eval('bsxfun(@times,transpose(modsigma),f)');


% Calculating the price index (P_is)
for i=1:s,
    origin_shifter(:,:,i)=eval('bsxfun(@times,secsize(:,i),pricingsigma(:,:,i))');
end
priceind=sum(origin_shifter,1);

exponentstore=1./minussigmat;
priceind_squeezed=squeeze(priceind);

for i=1:s,
    realpriceind(:,i)=priceind_squeezed(:,i).^(exponentstore(i));
end

% Calculating Trade shares
for i=1:s,
    tradeshare(:,:,i)=eval('bsxfun(@ldivide,priceind(:,:,i),origin_shifter(:,:,i))');
end
% Calculating Incomes
inc=wt.*L;
secexp=eval('bsxfun(@times,mu,inc)');
secinc=eval('bsxfun(@times,wt,Ls)');
% Calculating Market access

for i=1:s,
    tausigma(:,:,i)=eval('bsxfun(@power,tau(:,:,i,r),minussigmat(i))');
end

for i=1:s,
    realexp(:,i)=secexp(:,i)./priceind_squeezed(:,i);
    market_acc(:,i)=tausigma(:,:,i)*realexp(:,i);
end

 %% Checks on market access
% %Checking Market access by calculating from equilibrium condition
% for i=1:s,
%     divisor_1(:,i)=eval('bsxfun(@times,markupt(i),unitcost(:,i))');
%     divisor_2(:,i)=divisor_1(:,i).^minussigmat(i);
% end
% divisor=secsize.*divisor_2;
% market_acc_t=secinc./divisor;
% error_MA=market_acc_t-market_acc;
% 
% 
% 
% for i=1:s,
%     unitcostpower(:,i)=eval('bsxfun(@power,unitcost(:,i),minussigmat(i))');
% end
% 
% secinc(1,:)./secinc(2,:)-((unitcostpower(1,:)./unitcostpower(2,:)).*(secsize(1,:)./secsize(2,:)).*(market_acc(1,:)./market_acc(2,:)));
% 

