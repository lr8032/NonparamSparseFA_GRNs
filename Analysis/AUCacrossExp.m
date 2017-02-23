% plot the errorbars of mutiple samples' AUC values across experiments
%
%
clear all; close all; clc

% E(1)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter4000trial1_yeast1000.mat');
% E(2)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter4100trial1_yeast1000.mat');
% E(3)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3900trial1_yeast1000.mat');
% E(4)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3800trial1_yeast1000.mat');
% E(5)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3700trial1_yeast1000.mat');
% E(6)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3600trial1_yeast1000.mat');
% load RandTestYeast1000.mat

% E(1)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9200trial1_saureus.mat');
% E(2)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9300trial1_saureus.mat');
% E(3)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9400trial1_saureus.mat');
% E(4)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9500trial1_saureus.mat');
% E(5)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9600trial1_saureus.mat');
% E(6)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9700trial1_saureus.mat');
% E(7)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9800trial1_saureus.mat');
% E(8)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9900trial1_saureus.mat');
% E(9)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter10000trial1_saureus.mat');
% load /Users/rxlics/Desktop/Gene/RandTestSaureus1018.mat

E(1)=load('/Users/rxlics/Desktop/Gene/test_from_RC/test_ecoli1500/BNP_covreg_statsiter700trial3.mat');
E(2)=load('/Users/rxlics/Desktop/Gene/test_from_RC/test_ecoli1500/BNP_covreg_statsiter600trial3.mat');
E(3)=load('/Users/rxlics/Desktop/Gene/test_from_RC/test_ecoli1500/BNP_covreg_statsiter500trial3.mat');
% E(4)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3800trial1_yeast1000.mat');
% E(5)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3700trial1_yeast1000.mat');
% E(6)=load('/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter3600trial1_yeast1000.mat');
load RandTestEcoli1500.mat


% load C:\MATLAB\R2009a\work\Projects\Gene\Results\yeast1500\BNP_covreg_statsiter6500_yeast1500.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast1500.mat

% load C:\MATLAB\R2009a\work\Projects\Gene\Results\Saureus\BNP_covreg_statsiter3000trial1_Saureus.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestEcoli1500.mat
%%
[numG,b]=size(testGID);
% trueLabel=[ones(1,numG/4),zeros(1,numG/4)];
trueLabel=[ones(1,length(rowPos)),zeros(1,length(rowNeg))];
totAUC=[];
totSave=length(E);

for i=1:totSave
for saveIter=1:10;
    theta=E(i).Stats(saveIter).theta;
    zeta=E(i).Stats(saveIter).zeta;
    invSig_vec=E(i).Stats(saveIter).invSig_vec;
    Sigma=zeros(numG,numG);
    totExp=size(zeta,3);

    for tt=1:totExp
        Sigma(:,:)=theta*zeta(:,:,tt)*zeta(:,:,tt)'*theta' + diag(1./invSig_vec);
        for j=1:2:numG-1
            myScore((j+1)/2,tt)=(Sigma(j,j+1)+Sigma(j+1,j))/2;
        end
        [myFPrate(:,tt), myTPrate(:,tt), myAUC(saveIter,tt)]=plotROCcurve(myScore(:,tt),trueLabel);

    end
end
totAUC=[totAUC;myAUC];
end

errorbar(1:totExp, mean(totAUC),std(totAUC),'xb')
xlim([-1,totExp+2])

[B,IX] = sort(mean(totAUC));

for ii=1
    ind=IX(end-ii+1);
    Sigma=theta*zeta(:,:,ind)*zeta(:,:,ind)'*theta' + diag(1./invSig_vec);
    figure
    imagesc(Sigma)
    xlabel('Genes','fontsize',10,'fontweight','b') 
    ylabel('Genes','fontsize',10,'fontweight','b')
    title(['Covariance under condition: ' num2str(ind)],'fontsize',10,'fontweight','b')


%     for j=1:2:numG-1
%         indScore((j+1)/2,ind)=(Sigma(j,j+1)+Sigma(j+1,j))/2;
%     end
%     [indFPrate(:,ind), indTPrate(:,ind), indAUC(:,ind)]=plotROCcurve(indScore(:,ind),trueLabel);
%     figure
%     plot(myFPrate(:,ind),myTPrate(:,ind),'r','linewidth',2)
%     hold on
%     plot([0,1],[0,1],'k--')
%     axis([0,1,0,1])
%     title(['ROC under condition: ' num2str(ind)],'fontsize',10,'fontweight','b')
%     xlabel('False positive rate','fontsize',10,'fontweight','b')
%     ylabel('True positive rate','fontsize',10,'fontweight','b')
end