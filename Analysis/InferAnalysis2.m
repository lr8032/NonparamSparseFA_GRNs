
clear all; close all; clc

% load C:\MATLAB\R2009a\work\Projects\Gene\Results\BNP_covreg_statsiter10000trial1.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast50.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\yeast1000\BNP_covreg_statsiter4000trial1.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast1000.mat

% load C:\MATLAB\R2009a\work\Projects\Gene\Results\Saureus\BNP_covreg_statsiter9600trial1_saureus.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestSaureus1018.mat

load /Users/rxlics/Desktop/Gene/test_from_RC/test_ecoli1500/BNP_covreg_statsiter800trial3.mat
load RandTestEcoli1500.mat

%%
saveIter=1;
theta=Stats(saveIter).theta;
zeta=Stats(saveIter).zeta;
invSig_vec=Stats(saveIter).invSig_vec;
psi=Stats(saveIter).psi;
eta=Stats(saveIter).eta;

[numG,b]=size(theta);
[k,numExp]=size(psi);

Sigma=zeros(numG,numG);
sumSigma=zeros(numG,numG);
mu=zeros(numG,numExp);
totExp=size(zeta,3);
for tt=1:totExp
    mu(:,tt)=theta*zeta(:,:,tt)*eta(:,tt);
    Sigma(:,:)=theta*zeta(:,:,tt)*zeta(:,:,tt)'*theta' + diag(1./invSig_vec);
%     corr34(tt)=Sigma(17,18); % 5-6
%     tmp=Sigma(17,:);
%     tmp(17)=[];
%     tmp(18)=[];
%     corr1rest(tt)=mean(tmp);
    var(:,tt)=diag(Sigma);
end

mGExp=mean(testGExp);
stdGExp=std(testGExp);

plot_fill = @(x,lower,upper,color,alphaValue) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color,'FaceAlpha',alphaValue);

figure; hold on
% plot(testGExp','.b','markersize',1)
plot(mGExp,'b','linewidth',2)
plot_fill(1:totExp,mGExp-stdGExp,mGExp+stdGExp,'b',0.3)
plot(mean(mu),'r','linewidth',1.5)
plot_fill(1:totExp,mean(mu)-mean(sqrt(var)),mean(mu)+mean(sqrt(var)),'r',0.3)
xlabel('Experiment Conditions','fontsize',12)
ylabel('Expression Level','fontsize',12)
xlim([1 totExp])
% legend('Empirical Stats',' ','Model Est',' ')

% figure
% plot(corr34,'r')
% hold on
% plot(corr1rest,'b')

% figure; hold on
% boxplot(testGExp','plotstyle','compact')
% plot(median(testGExp'),'b','linewidth',2)
% plot(mean(mu),'r','linewidth',1.5)
% plot_fill(1:536,mean(mu)-mean(sqrt(var)),mean(mu)+mean(sqrt(var)),'r',0.3)
% xlim([1 536])

% figure; hold on
% plot(mu','r','linewidth',1.5)
% plot(median(testGExp),'b','linewidth',2)







