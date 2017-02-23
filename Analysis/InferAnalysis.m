function InferAnalysis()
clear all; close all; clc

% load /Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter5100trial1.mat
% load /Users/rxlics/Desktop/Gene/RandTestYeast1000.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\yeast1500\BNP_covreg_statsiter6500_yeast1500.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast1500.mat

% load /Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter1000trial1.mat
% load /Users/rxlics/Desktop/Gene/ecoli_from_Mac_Pro/BNP_covreg_statsiter1000trial1.mat
% load /Users/rxlics/Desktop/Gene/test_from_RC/test_ecoli1500/BNP_covreg_statsiter800trial3.mat
% load /Users/rxlics/Desktop/Gene/RandTestEcoli1500.mat

load /Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test/BNP_covreg_statsiter9600trial1_saureus.mat
load /Users/rxlics/Desktop/Gene/RandTestSaureus1018.mat

% ciA=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Community integration/DREAM5_NetworkInference_Community_Network2.txt');
% meta1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Meta1/DREAM5_NetworkInference_Meta1_Network2.txt');
% meta2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Meta2/DREAM5_NetworkInference_Meta2_Network2.txt');
% meta3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Meta3/DREAM5_NetworkInference_Meta3_Network2.txt');
% meta4A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Meta4/DREAM5_NetworkInference_Meta4_Network2.txt');
% meta5A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Meta5/DREAM5_NetworkInference_Meta5_Network2.txt');
% corr1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Correlation1/DREAM5_NetworkInference_Correlation1_Network2.txt');
% corr2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/Correlation2/DREAM5_NetworkInference_Correlation2_Network2.txt');
% corr3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/Correlation3/DREAM5_NetworkInference_Correlation3_Network2.txt');
% mi1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/MI1/DREAM5_NetworkInference_MI1_Network2.txt');
% mi2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/MI2/DREAM5_NetworkInference_MI2_Network2.txt');
% mi3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/MI3/DREAM5_NetworkInference_MI3_Network2.txt');
% mi4A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/MI4/DREAM5_NetworkInference_MI4_Network2.txt');
% mi5A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/MI5/DREAM5_NetworkInference_MI5_Network2.txt');
% reg1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression1/DREAM5_NetworkInference_Regression1_Network2.txt');
% reg2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression2/DREAM5_NetworkInference_Regression2_Network2.txt');
% reg3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression3/DREAM5_NetworkInference_Regression3_Network2.txt');
% reg4A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression4/DREAM5_NetworkInference_Regression4_Network2.txt');
% reg5A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression5/DREAM5_NetworkInference_Regression5_Network2.txt');
% reg6A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression6/DREAM5_NetworkInference_Regression6_Network2.txt');
% reg7A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Regression7/DREAM5_NetworkInference_Regression7_Network2.txt');
% reg8A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Off-the-shelf methods/Regression8/DREAM5_NetworkInference_Regression8_Network2.txt');
% bay1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian1/DREAM5_NetworkInference_Bayesian1_Network2.txt');
% bay2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian2/DREAM5_NetworkInference_Bayesian2_Network2.txt');
% bay3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian3/DREAM5_NetworkInference_Bayesian3_Network2.txt');
% bay4A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian4/DREAM5_NetworkInference_Bayesian4_Network2.txt');
% bay5A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian5/DREAM5_NetworkInference_Bayesian5_Network2.txt');
% bay6A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Bayesian6/DREAM5_NetworkInference_Bayesian6_Network2.txt');
% other1A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other1/DREAM5_NetworkInference_Other1_Network2.txt');
% other2A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other2/DREAM5_NetworkInference_Other2_Network2.txt');
% other3A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other3/DREAM5_NetworkInference_Other3_Network2.txt');
% other4A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other4/DREAM5_NetworkInference_Other4_Network2.txt');
% other5A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other5/DREAM5_NetworkInference_Other5_Network2.txt');
% other6A=importdata('/Users/rxlics/Desktop/Gene/DREAM5_methods/Network_predictions/Challenge participants/Other6/DREAM5_NetworkInference_Other6_Network2.txt');


%%
[numG,b]=size(testGID);
% trueLabel=[ones(1,numG/4),zeros(1,numG/4)];
trueLabel=[ones(1,length(rowPos)),zeros(1,length(rowNeg))];

saveIter=1;
theta=Stats(saveIter).theta;
zeta=Stats(saveIter).zeta;
invSig_vec=Stats(saveIter).invSig_vec;
Sigma=zeros(numG,numG);
sumSigma=zeros(numG,numG);

totExp=size(zeta,3);
% figure
for tt=1:totExp
%     subplot(1,3,1)
%     imagesc(theta)
%     subplot(1,3,2)
%     imagesc(zeta(:,:,tt))
%     subplot(1,3,3)
%     imagesc(diag(1./invSig_vec))
%     pause
    
    Sigma(:,:)=theta*zeta(:,:,tt)*zeta(:,:,tt)'*theta' + diag(1./invSig_vec);
%     subplot(1,2,1)
%     imagesc(Sigma(:,:))
%     xlabel('Genes','fontsize',10,'fontweight','b') 
%     ylabel('Genes','fontsize',10,'fontweight','b')
%     title(['Covariance under condition: ' num2str(tt)],'fontsize',10,'fontweight','b')
    
    for j=1:2:numG-1
        myScore((j+1)/2,tt)=(Sigma(j,j+1)+Sigma(j+1,j))/2;
    end
    [myFPrate(:,tt), myTPrate(:,tt), myAUC(:,tt)]=plotROCcurve(myScore(:,tt),trueLabel);
%     subplot(1,2,2)
%     plot(myFPrate(:,tt),myTPrate(:,tt),'r','linewidth',2)
%     hold on
%     plot([0,1],[0,1],'k--')
%     axis([0,1,0,1])
%     title(['ROC under condition: ' num2str(tt)],'fontsize',10,'fontweight','b')
%     xlabel('False positive rate','fontsize',10,'fontweight','b')
%     ylabel('True positive rate','fontsize',10,'fontweight','b')
%     if tt==1
%         pause
%     end
%     pause
%     hold off
%     sumSigma=sumSigma+Sigma;
end


[B,IX] = sort(myAUC);

for ii=1
    ind=IX(end-ii+1);
    Sigma=theta*zeta(:,:,ind)*zeta(:,:,ind)'*theta' + diag(1./invSig_vec);
    figure
    imagesc(Sigma)
    xlabel('Genes','fontsize',10,'fontweight','b') 
    ylabel('Genes','fontsize',10,'fontweight','b')
    title(['Covariance under condition: ' num2str(ind)],'fontsize',10,'fontweight','b')


    for j=1:2:numG-1
        indScore((j+1)/2,ind)=(Sigma(j,j+1)+Sigma(j+1,j))/2;
    end
    [indFPrate(:,ind), indTPrate(:,ind), indAUC]=plotROCcurve(indScore(:,ind),trueLabel);
    save DREAM5_AUC.mat indAUC -append
    figure
    plot(indFPrate(:,ind),indTPrate(:,ind),'r','linewidth',3)
    hold on
    plot([0,1],[0,1],'k--','linewidth',2)
    axis([0,1,0,1])
    title(['ROC under condition: ' num2str(ind)],'fontsize',10,'fontweight','b')
    xlabel('False positive rate','fontsize',10,'fontweight','b')
    ylabel('True positive rate','fontsize',10,'fontweight','b')
end

%%
% pertubAUC=[myAUC(:,119:120),myAUC(:,125:148),myAUC(:,153:177),myAUC(:,180:203),myAUC(:,206:221),myAUC(:,225:227),myAUC(:,231:233),myAUC(:,236:241),myAUC(:,245:247),myAUC(:,251:255),myAUC(:,262:271),myAUC(:,282:291),myAUC(:,303:361),myAUC(:,369:378),myAUC(:,496:536)];
% unpertubAUC=[myAUC(:,1:118),myAUC(:,121:124),myAUC(:,149:152),myAUC(:,178:179),myAUC(:,204:205),myAUC(:,222:224),myAUC(:,228:230),myAUC(:,234:235),myAUC(:,242:244),myAUC(:,248:250),myAUC(:,256:261),myAUC(:,272:281),myAUC(:,292:302),myAUC(:,362:368),myAUC(:,379:495)];
% figure
% hold on
% pertubH=bar(1,mean(pertubAUC));
% set(pertubH,'FaceColor','r');
% unpertubH=bar(2,mean(unpertubAUC));
% set(unpertubH,'FaceColor','b');
% errorbar(1:2,[mean(pertubAUC),mean(unpertubAUC)],[std(pertubAUC),std(unpertubAUC)],'k.','linewidth',2)
% ax=gca;
% XTick={' ' 'Pertubed' 'Unpertubed' ' '};
% set(ax,'XLim',[0,3])
% set(ax,'XTick',[0:3])
% set(ax,'XTickLabel',XTick)
% hold off
% 
% deleteAUC=[myAUC(:,382:384),myAUC(:,388:396),myAUC(:,399:411),myAUC(:,418:435),myAUC(:,466:495),myAUC(:,504:536)];
% undeleteAUC=[myAUC(:,1:381),myAUC(:,385:387),myAUC(:,397:398),myAUC(:,412:417),myAUC(:,436:465),myAUC(:,496:503)];
% figure
% hold on
% deleteH=bar(1,mean(deleteAUC));
% set(deleteH,'FaceColor','r');
% undeleteH=bar(2,mean(undeleteAUC));
% set(undeleteH,'FaceColor','b');
% errorbar(1:2,[mean(deleteAUC),mean(undeleteAUC)],[std(deleteAUC),std(undeleteAUC)],'k.','linewidth',2)
% ax=gca;
% XTick={' ' 'Deletion' 'Non-deletion' ' '};
% set(ax,'XLim',[0,3])
% set(ax,'XTick',[0:3])
% set(ax,'XTickLabel',XTick)
% hold off


%%
[ciScore,ciLabel,ciPair]=DreamROCofTest(ciA,testGID);
[ciFPrate, ciTPrate, ciAUC]=plotROCcurve(ciScore,ciLabel);
ciH=plot(ciFPrate,ciTPrate,'color',[0.6 0.6 0.6],'linewidth',1.5);
save DREAM5_AUC.mat ciAUC -append



[meta1Score,meta1Label,meta1Pair]=DreamROCofTest(meta1A,testGID);
[meta1FPrate, meta1TPrate, meta1AUC]=plotROCcurve(meta1Score,meta1Label);
meta1H=plot(meta1FPrate,meta1TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat meta1AUC -append
[meta2Score,meta2Label,meta2Pair]=DreamROCofTest(meta2A,testGID);
[meta2FPrate, meta2TPrate, meta2AUC]=plotROCcurve(meta2Score,meta2Label);
meta2H=plot(meta2FPrate,meta2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat meta2AUC -append
[meta3Score,meta3Label,meta3Pair]=DreamROCofTest(meta3A,testGID);
[meta3FPrate, meta3TPrate, meta3AUC]=plotROCcurve(meta3Score,meta3Label);
meta3H=plot(meta3FPrate,meta3TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat meta3AUC -append
[meta4Score,meta4Label,meta4Pair]=DreamROCofTest(meta4A,testGID);
[meta4FPrate, meta4TPrate, meta4AUC]=plotROCcurve(meta4Score,meta4Label);
meta4H=plot(meta4FPrate,meta4TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat meta4AUC -append
[meta5Score,meta5Label,meta5Pair]=DreamROCofTest(meta5A,testGID);
[meta5FPrate, meta5TPrate, meta5AUC]=plotROCcurve(meta5Score,meta5Label);
meta5H=plot(meta5FPrate,meta5TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat meta5AUC -append


[reg1Score,reg1Label,reg1Pair]=DreamROCofTest(reg1A,testGID);
[reg1FPrate, reg1TPrate, reg1AUC]=plotROCcurve(reg1Score,reg1Label);
reg1H=plot(reg1FPrate,reg1TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg1AUC -append
[reg2Score,reg2Label,reg2Pair]=DreamROCofTest(reg2A,testGID);
[reg2FPrate, reg2TPrate, reg2AUC]=plotROCcurve(reg2Score,reg2Label);
reg2H=plot(reg2FPrate,reg2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg2AUC -append
[reg3Score,reg3Label,reg3Pair]=DreamROCofTest(reg3A,testGID);
[reg3FPrate, reg3TPrate, reg3AUC]=plotROCcurve(reg3Score,reg3Label);
reg3H=plot(reg3FPrate,reg3TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg3AUC -append
[reg4Score,reg4Label,reg4Pair]=DreamROCofTest(reg4A,testGID);
[reg4FPrate, reg4TPrate, reg4AUC]=plotROCcurve(reg4Score,reg4Label);
reg4H=plot(reg4FPrate,reg4TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg4AUC -append
[reg5Score,reg5Label,reg5Pair]=DreamROCofTest(reg5A,testGID);
[reg5FPrate, reg5TPrate, reg5AUC]=plotROCcurve(reg5Score,reg5Label);
reg5H=plot(reg5FPrate,reg5TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg5AUC -append
[reg6Score,reg6Label,reg6Pair]=DreamROCofTest(reg6A,testGID);
[reg6FPrate, reg6TPrate, reg6AUC]=plotROCcurve(reg6Score,reg6Label);
reg6H=plot(reg6FPrate,reg6TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg6AUC -append
[reg7Score,reg7Label,reg7Pair]=DreamROCofTest(reg7A,testGID);
[reg7FPrate, reg7TPrate, reg7AUC]=plotROCcurve(reg7Score,reg7Label);
reg7H=plot(reg7FPrate,reg7TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg7AUC -append
[reg8Score,reg8Label,reg8Pair]=DreamROCofTest(reg8A,testGID);
[reg8FPrate, reg8TPrate, reg8AUC]=plotROCcurve(reg8Score,reg8Label);
reg8H=plot(reg8FPrate,reg8TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat reg8AUC -append



% [corr1Score,corr1Label,corr1Pair]=DreamROCofTest(corr1A,testGID);
% [corr1FPrate, corr1TPrate, corr1AUC]=plotROCcurve(corr1Score,corr1Label);
% corr1H=plot(corr1FPrate,corr1TPrate,'color',[.6 .6 .6],'linewidth',1.5);
% save DREAM5_AUC.mat corr1AUC -append
[corr2Score,corr2Label,corr2Pair]=DreamROCofTest(corr2A,testGID);
[corr2FPrate, corr2TPrate, corr2AUC]=plotROCcurve(corr2Score,corr2Label);
corr2H=plot(corr2FPrate,corr2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat corr2AUC -append
[corr3Score,corr3Label,corr3Pair]=DreamROCofTest(corr3A,testGID);
[corr3FPrate, corr3TPrate, corr3AUC]=plotROCcurve(corr3Score,corr3Label);
corr3H=plot(corr3FPrate,corr3TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat corr3AUC -append



[mi1Score,mi1Label,mi1Pair]=DreamROCofTest(mi1A,testGID);
[mi1FPrate, mi1TPrate, mi1AUC]=plotROCcurve(mi1Score,mi1Label);
mi1H=plot(mi1FPrate,mi1TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat mi1AUC -append
[mi2Score,mi2Label,mi2Pair]=DreamROCofTest(mi2A,testGID);
[mi2FPrate, mi2TPrate, mi2AUC]=plotROCcurve(mi2Score,mi2Label);
mi2H=plot(mi2FPrate,mi2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat mi2AUC -append
% [mi3Score,mi3Label,mi3Pair]=DreamROCofTest(mi3A,testGID);
% [mi3FPrate, mi3TPrate, mi3AUC]=plotROCcurve(mi3Score,mi3Label);
% mi3H=plot(mi3FPrate,mi3TPrate,'color',[.6 .6 .6],'linewidth',1.5);
% save DREAM5_AUC.mat mi3AUC -append
% [mi4Score,mi4Label,mi4Pair]=DreamROCofTest(mi4A,testGID);
% [mi4FPrate, mi4TPrate, mi4AUC]=plotROCcurve(mi4Score,mi4Label);
% mi4H=plot(mi4FPrate,mi4TPrate,'color',[.6 .6 .6],'linewidth',1.5);
% save DREAM5_AUC.mat mi4AUC -append
[mi5Score,mi5Label,mi5Pair]=DreamROCofTest(mi5A,testGID);
[mi5FPrate, mi5TPrate, mi5AUC]=plotROCcurve(mi5Score,mi5Label);
mi5H=plot(mi5FPrate,mi5TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat mi5AUC -append
% figure
% bar([mean(myAUC),ciAUC,corr2AUC,corr3AUC])
% mydata=rand(1,5);
% figure(1)
% hold on
% for i = 1:length(mydata)
%     h=bar(i,mydata(i));
%     if mydata(i) < 0.2
%         set(h,'FaceColor','k');
%     elseif mydata(i) < 0.6
%         set(h,'FaceColor','b');
%     else
%         set(h,'FaceColor','r');
%     end
% end
% hold off




% [bay1Score,bay1Label,bay1Pair]=DreamROCofTest(bay1A,testGID);
% [bay1FPrate, bay1TPrate, bay1AUC]=plotROCcurve(bay1Score,bay1Label);
% bay1H=plot(bay1FPrate,bay1TPrate,'color',[.6 .6 .6],'linewidth',2);
% save DREAM5_AUC.mat bay1AUC -append
[bay2Score,bay2Label,bay2Pair]=DreamROCofTest(bay2A,testGID);
[bay2FPrate, bay2TPrate, bay2AUC]=plotROCcurve(bay2Score,bay2Label);
bay2H=plot(bay2FPrate,bay2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat bay2AUC -append
% [bay3Score,bay3Label,bay3Pair]=DreamROCofTest(bay3A,testGID);
% [bay3FPrate, bay3TPrate, bay3AUC]=plotROCcurve(bay3Score,bay3Label);
% bay3H=plot(bay3FPrate,bay3TPrate,'color',[.6 .6 .6],'linewidth',2);
% save DREAM5_AUC.mat bay3AUC -append
[bay4Score,bay4Label,bay4Pair]=DreamROCofTest(bay4A,testGID);
[bay4FPrate, bay4TPrate, bay4AUC]=plotROCcurve(bay4Score,bay4Label);
bay4H=plot(bay4FPrate,bay4TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat bay4AUC -append
[bay5Score,bay5Label,bay5Pair]=DreamROCofTest(bay5A,testGID);
[bay5FPrate, bay5TPrate, bay5AUC]=plotROCcurve(bay5Score,bay5Label);
bay5H=plot(bay5FPrate,bay5TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat bay5AUC -append
% [bay6Score,bay6Label,bay6Pair]=DreamROCofTest(bay6A,testGID);
% [bay6FPrate, bay6TPrate, bay6AUC]=plotROCcurve(bay6Score,bay6Label);
% bay6H=plot(bay6FPrate,bay6TPrate,'color',[.6 .6 .6],'linewidth',1.5);
% save DREAM5_AUC.mat bay6AUC -append


[other1Score,other1Label,other1Pair]=DreamROCofTest(other1A,testGID);
[other1FPrate, other1TPrate, other1AUC]=plotROCcurve(other1Score,other1Label);
other1H=plot(other1FPrate,other1TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat other1AUC -append
[other2Score,other2Label,other2Pair]=DreamROCofTest(other2A,testGID);
[other2FPrate, other2TPrate, other2AUC]=plotROCcurve(other2Score,other2Label);
other2H=plot(other2FPrate,other2TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat other2AUC -append
[other3Score,other3Label,other3Pair]=DreamROCofTest(other3A,testGID);
[other3FPrate, other3TPrate, other3AUC]=plotROCcurve(other3Score,other3Label);
other3H=plot(other3FPrate,other3TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat other3AUC -append
% [other4Score,other4Label,other4Pair]=DreamROCofTest(other4A,testGID);
% [other4FPrate, other4TPrate, other4AUC]=plotROCcurve(other4Score,other4Label);
% other4H=plot(other4FPrate,other4TPrate,'color',[.6 .6 .6],'linewidth',1.5);
% save DREAM5_AUC.mat other4AUC -append
[other5Score,other5Label,other5Pair]=DreamROCofTest(other5A,testGID);
[other5FPrate, other5TPrate, other5AUC]=plotROCcurve(other5Score,other5Label);
other5H=plot(other5FPrate,other5TPrate,'color',[.6 .6 .6],'linewidth',1.5);
save DREAM5_AUC.mat other5AUC -append
% [other6Score,other6Label,other6Pair]=DreamROCofTest(other6A,testGID);
% [other6FPrate, other6TPrate, other6AUC]=plotROCcurve(other6Score,other6Label);
% other6H=plot(other6FPrate,other6TPrate,'color',[.6 .6 .6],'linewidth',2);
% save DREAM5_AUC.mat other6AUC -append



plot(indFPrate(:,ind),indTPrate(:,ind),'r','linewidth',3)
reg5H=plot(reg5FPrate,reg5TPrate,'color',[.1 .1 .1],'linewidth',3);
plot([0,1],[0,1],'k--','linewidth',2)
    
    
% for i=1:2:a-1
%    score((i+1)/2)=(abs(rho(i,i+1))+abs(rho(i+1,i)))/2; 
% end
% plotROCcurve(score,trueLabel,1)

% ciScore=[];
% ciLabel=[];
% ciPair=[];
% for j=1:length(ciA.data)
%     ind_tf=strmatch(ciA.textdata{j,1},testGID,'exact');
%     if ind_tf~=0 % find a match gene, which may be repeated in testGID
%         for k=1:length(ind_tf) % for each repetition of the gene
%             tf_len=length(ciA.textdata{j,2});
%             tg_len=length(strtrim(testGID(ind_tf(k)-1+2*mod(ind_tf(k),2),:)));
%             max_len=max([tf_len,tg_len]);
%             if strncmp(ciA.textdata{j,2},testGID(ind_tf(k)-1+2*mod(ind_tf(k),2),:),max_len) % find its pair also matches
%                 ciScore=[ciScore,ciA.data(j)];
%                 if ind_tf(k)<=numG/2 % if the pair is in the first half, it's true
%                     ciLabel=[ciLabel,1];
%                 else
%                     ciLabel=[ciLabel,0];
%                 end
%                 if mod(ind_tf(k),2)==1
%                     ciPair=[ciPair,ind_tf(k)];
%                 else
%                     ciPair=[ciPair,ind_tf(k)+1];
%                 end
%             end
%         end
%     end
% end





% for i=1:536
%     for j=1:2:numG-1
%         myScore((j+1)/2,i)=(Sigma(j,j+1,i)+Sigma(j+1,j,i))/2;
%     end
%     [myFPrate(:,i), myTPrate(:,i), myAUC(:,i)]=plotROCcurve(myScore(:,i),trueLabel);
%     plot(myFPrate(:,i),myTPrate(:,i),'r','linewidth',2)
%     hold on
%     plot(ciFPrate,ciTPrate,'b','linewidth',2)
%     plot([0,1],[0,1],'k--')
%     axis([0,1,0,1])
%     title(['ROC under condition: ' num2str(i)],'FontSize',15)
%     xlabel('False positive rate','fontsize',10,'fontweight','b')
%     ylabel('True positive rate','fontsize',10,'fontweight','b')
%     pause
%     hold off
% end
%   myUniFPrate=unique(myFPrate);
% for t=1:length(myUniFPrate)
%    myMTPrate(t)=mode(myTPrate(myFPrate==myUniFPrate(t))); 
% end
% plot(myUniFPrate,myMTPrate,'r','linewidth',2)
% averS=sum(myScore,2)/536;
% plotROCcurve(averS,trueLabel,1,'r')
% [X,Y,THRE,AUC,OPTROCPT,SUBY,SUBYNAMES]=perfcurve(trueLabel',averS,1);
% plot(X,Y)
% hold on
% plot(mean(myFPrate,2),mean(myTPrate,2),'r','linewidth',2)
% plot(mean(myFPrate,2),mean(myTPrate,2)-std(myTPrate,1,2),'--','color',[0.8 0 0],'linewidth',2)
% plot(mean(myFPrate,2),mean(myTPrate,2)+std(myTPrate,1,2),'--','color',[0.8 0 0],'linewidth',2)
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')
% grid on
end






%%
function [Score,Label,Pair]=DreamROCofTest(A,testGID)
[numG,b]=size(testGID);
Score=[];
Label=[];
Pair=[];
for j=1:length(A.data)
    ind_tf=strmatch(A.textdata{j,1},testGID,'exact');
    if ind_tf~=0 % find a match gene, which may be repeated in testGID
        for k=1:length(ind_tf) % for each repetition of the gene
            tf_len=length(A.textdata{j,2});
            tg_len=length(strtrim(testGID(ind_tf(k)-1+2*mod(ind_tf(k),2),:))); % if ind_tf is even, its partner is ind_tf-1; if ind_tf is odd, its partner is ind_tf+1
            max_len=max([tf_len,tg_len]);
            if strncmp(A.textdata{j,2},testGID(ind_tf(k)-1+2*mod(ind_tf(k),2),:),max_len) % find its pair also matches
                Score=[Score,A.data(j)];
                if ind_tf(k)<=numG/2 % if the pair is in the first half, it's true
                    Label=[Label,1];
                else
                    Label=[Label,0];
                end
                if mod(ind_tf(k),2)==1
                    Pair=[Pair,ind_tf(k)];
                else
                    Pair=[Pair,ind_tf(k)+1];
                end
            end
        end
    end
end


end
