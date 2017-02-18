function InferAnalysis()
clear all; close all; clc

load C:\MATLAB\R2009a\work\Projects\Gene\Results\yeast1000\BNP_covreg_statsiter4600trial1.mat
load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast1000.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\yeast1500\BNP_covreg_statsiter6500_yeast1500.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestYeast1500.mat

% load C:\MATLAB\R2009a\work\Projects\Gene\Results\Saureus\BNP_covreg_statsiter7500trial1.mat
% load C:\MATLAB\R2009a\work\Projects\Gene\Results\RandTestSaureus1018.mat

% ciA=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Community integration\DREAM5_NetworkInference_Community_Network2.txt');
% corr1A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Correlation1\DREAM5_NetworkInference_Correlation1_Network2.txt');
% corr2A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Off-the-shelf methods\Correlation2\DREAM5_NetworkInference_Correlation2_Network2.txt');
% corr3A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Off-the-shelf methods\Correlation3\DREAM5_NetworkInference_Correlation3_Network2.txt');
% mi1A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Off-the-shelf methods\MI1\DREAM5_NetworkInference_MI1_Network2.txt');
% reg1A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Regression1\DREAM5_NetworkInference_Regression1_Network2.txt');
% meta1A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Meta1\DREAM5_NetworkInference_Meta1_Network2.txt');
% meta2A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Meta2\DREAM5_NetworkInference_Meta2_Network2.txt');
% meta3A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Meta3\DREAM5_NetworkInference_Meta3_Network2.txt');
% meta4A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Meta4\DREAM5_NetworkInference_Meta4_Network2.txt');
% meta5A=importdata('C:\MATLAB\R2009a\work\Projects\Gene\DREAM5_methods\Network_predictions\Challenge participants\Meta5\DREAM5_NetworkInference_Meta5_Network2.txt');

%%
[numG,b]=size(testGID);
trueLabel=[ones(1,numG/4),zeros(1,numG/4)];

saveIter=10;
theta=Stats(saveIter).theta;
zeta=Stats(saveIter).zeta;
invSig_vec=Stats(saveIter).invSig_vec;
Sigma=zeros(numG,numG);
% sumSigma=zeros(numG,numG);

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
    hold off
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
    [indFPrate(:,ind), indTPrate(:,ind), indAUC(:,ind)]=plotROCcurve(indScore(:,ind),trueLabel);
    figure
    plot(myFPrate(:,ind),myTPrate(:,ind),'r','linewidth',2)
    hold on
    plot([0,1],[0,1],'k--')
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
% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[meta1Score,meta1Label,meta1Pair]=DreamROCofTest(meta1A,testGID);
[meta1FPrate, meta1TPrate, meta1AUC]=plotROCcurve(meta1Score,meta1Label);
meta1H=plot(meta1FPrate,meta1TPrate,'b','linewidth',2);

[meta2Score,meta2Label,meta2Pair]=DreamROCofTest(meta2A,testGID);
[meta2FPrate, meta2TPrate, meta2AUC]=plotROCcurve(meta2Score,meta2Label);
meta2H=plot(meta2FPrate,meta2TPrate,'b','linewidth',2);

[meta3Score,meta3Label,meta3Pair]=DreamROCofTest(meta3A,testGID);
[meta3FPrate, meta3TPrate, meta3AUC]=plotROCcurve(meta3Score,meta3Label);
meta3H=plot(meta3FPrate,meta3TPrate,'b','linewidth',2);

[meta4Score,meta4Label,meta4Pair]=DreamROCofTest(meta4A,testGID);
[meta4FPrate, meta4TPrate, meta4AUC]=plotROCcurve(meta4Score,meta4Label);
meta4H=plot(meta4FPrate,meta4TPrate,'b','linewidth',2);

[meta5Score,meta5Label,meta5Pair]=DreamROCofTest(meta5A,testGID);
[meta5FPrate, meta5TPrate, meta5AUC]=plotROCcurve(meta5Score,meta5Label);
meta5H=plot(meta5FPrate,meta5TPrate,'b','linewidth',2);



[ciScore,ciLabel,ciPair]=DreamROCofTest(ciA,testGID);
[ciFPrate, ciTPrate, ciAUC]=plotROCcurve(ciScore,ciLabel);
ciH=plot(ciFPrate,ciTPrate,'b','linewidth',2);

% sumSigma=sumSigma/totExp;
% 
% for i=1:length(ciPair)
%    myCoScore(i)=(sumSigma(ciPair(i),ciPair(i)+1)+sumSigma(ciPair(i)+1,ciPair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,ciLabel);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([ciH myCoH],{'community integration','my method'},4)



% 
% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[reg1Score,reg1Label,reg1Pair]=DreamROCofTest(reg1A,testGID);
[reg1FPrate, reg1TPrate, reg1AUC]=plotROCcurve(reg1Score,reg1Label);
reg1H=plot(reg1FPrate,reg1TPrate,'g','linewidth',2);

% for i=1:length(reg1Pair)
%    myCoScore(i)=(sumSigma(reg1Pair(i),reg1Pair(i)+1)+sumSigma(reg1Pair(i)+1,reg1Pair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,reg1Label);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([reg1H myCoH],{'regression1','my method'},4)




% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[corr1Score,corr1Label,corr1Pair]=DreamROCofTest(corr1A,testGID);
[corr1FPrate, corr1TPrate, corr1AUC]=plotROCcurve(corr1Score,corr1Label);
corr1H=plot(corr1FPrate,corr1TPrate,'g','linewidth',2);

% for i=1:length(corr1Pair)
%    myCoScore(i)=(sumSigma(corr1Pair(i),corr1Pair(i)+1)+sumSigma(corr1Pair(i)+1,corr1Pair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,corr1Label);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([corr1H myCoH],{'correlation1','my method'},4)



% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[corr2Score,corr2Label,corr2Pair]=DreamROCofTest(corr2A,testGID);
[corr2FPrate, corr2TPrate, corr2AUC]=plotROCcurve(corr2Score,corr2Label);
corr2H=plot(corr2FPrate,corr2TPrate,'g','linewidth',2);

% for i=1:length(corr2Pair)
%    myCoScore(i)=(sumSigma(corr2Pair(i),corr2Pair(i)+1)+sumSigma(corr2Pair(i)+1,corr2Pair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,corr2Label);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([corr2H myCoH],{'correlation2','my method'},4)



% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[corr3Score,corr3Label,corr3Pair]=DreamROCofTest(corr3A,testGID);
[corr3FPrate, corr3TPrate, corr3AUC]=plotROCcurve(corr3Score,corr3Label);
corr3H=plot(corr3FPrate,corr3TPrate,'g','linewidth',2);

% for i=1:length(corr3Pair)
%    myCoScore(i)=(sumSigma(corr3Pair(i),corr3Pair(i)+1)+sumSigma(corr3Pair(i)+1,corr3Pair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,corr3Label);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([corr3H myCoH],{'correlation3','my method'},4)


% figure
% hold on
% plot([0,1],[0,1],'k--')
% axis([0,1,0,1])
% xlabel('False positive rate','fontsize',10,'fontweight','b')
% ylabel('True positive rate','fontsize',10,'fontweight','b')

[mi1Score,mi1Label,mi1Pair]=DreamROCofTest(mi1A,testGID);
[mi1FPrate, mi1TPrate, mi1AUC]=plotROCcurve(mi1Score,mi1Label);
mi1H=plot(mi1FPrate,mi1TPrate,'g','linewidth',2);

% for i=1:length(mi1Pair)
%    myCoScore(i)=(sumSigma(mi1Pair(i),mi1Pair(i)+1)+sumSigma(mi1Pair(i)+1,mi1Pair(i)))/2; 
% end
% [myCoFPrate, myCoTPrate, myCoAUC]=plotROCcurve(myCoScore,mi1Label);
% myCoH=plot(myCoFPrate,myCoTPrate,'r','linewidth',2);
% legend([mi1H myCoH],{'mi3','my method'},4)



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
