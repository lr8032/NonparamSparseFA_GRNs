% plot the gene expressions and the regulatory network structures.
%
% written by RL 01/17

clear all; close all; clc

% load '/Users/rxlics/Desktop/Gene/BNPCovReg_toolbox/test on 100 Gene Pairs/BNP_covreg_statsiter10000trial1.mat'
% load '/Users/rxlics/Desktop/Gene/RandTestYeast50.mat'
% saveIter=1;
% [numG,b]=size(Stats(saveIter).phi);
% 
% theta=Stats(saveIter).theta;
% zeta=Stats(saveIter).zeta;
% invSig_vec=Stats(saveIter).invSig_vec;
% Sigma=zeros(numG,numG);
% 
% tt=300; %1:536
% Sigma(:,:)=theta*zeta(:,:,tt)*zeta(:,:,tt)'*theta' + diag(1./invSig_vec);
% 
% % generate the graph
% numPlot=6;
% % names={'YOR335C' 'YBR049C' 'YOL040C' 'YPR104C' 'YNL063W' 'YJR060W'};% 'YPR013C' 'YMR016C'};
% G=graph(Sigma(1:numPlot,1:numPlot).'| Sigma(1:numPlot,1:numPlot),'omitselfloops');
% 
% 
% % find weights on the connections
% weights=[];
% for i=1:numPlot-1
%     weights=[weights Sigma(i,i+1:numPlot)];
% end
% weights(1,weights<0)=-weights(1,weights<0);
% 
% % eliminate some edges
% thr=mean(weights);
% G=rmedge(G,find(weights<thr));
% weights(weights<thr)=[];
% 
% 
% figure
% colormap(gray)
% h=plot(G,'edgecdata',weights,'linewidth',3,'markersize',10,'layout','circle');
% highlight(h,1,'NodeColor','r')
% highlight(h,2,'NodeColor','g')
% highlight(h,3,'NodeColor','b')
% highlight(h,4,'NodeColor','m')
% highlight(h,5,'NodeColor','y')
% highlight(h,6,'NodeColor','c')
% hold on
% colormap(parula)
% plot(G,'edgecdata',zeros(1,numPlot),'linewidth',1,'markersize',6,'layout','circle','nodecdata',1:numPlot)

% figure
% plot(testGExp(1,:)','.','markersize',8,'color','r')
% hold on
% plot(testGExp(2,:)','.','markersize',8,'color','g')
% plot(testGExp(3,:)','.','markersize',8,'color','b')
% plot(testGExp(4,:)','.','markersize',8,'color','m')
% plot(testGExp(5,:)','.','markersize',8,'color','y')
% plot(testGExp(6,:)','.','markersize',8,'color','c')
% 
% xlim([0,538])
% ylim([min(min(testGExp(1:numPlot,:)))-0.5,max(max(testGExp(1:numPlot,:)))+0.5])



%%
% load EstVsStand_test_yeast1000.mat
% load EstVsStand_test_saureus.mat
load EstVsStand_test_ecoli07.mat

% load nonRepGoldStand_yeast.mat
% load nonRepGoldStand_test_saureus.mat
% load nonRepGoldStand_test_ecoli.mat

GG=graph(str2num(cGeneA(1:end,2:end)),str2num(cGeneB(1:end,2:end)));

% h=plot(GG,'layout','force','iteration',600,'nodecolor','k');
h=plot(GG,'layout','force','iteration',600,'nodecolor','k','edgecolor','b','linewidth',1.5);

% fnG=graph(str2num(cGeneA(cIndicator==1,2:end)),str2num(cGeneB(cIndicator==1,2:end)));
% tpG=graph(str2num(cGeneA(cIndicator==4,2:end)),str2num(cGeneB(cIndicator==4,2:end))); 
% fpG=graph(str2num(cGeneA(cIndicator==3,2:end)),str2num(cGeneB(cIndicator==3,2:end)));
% % highlight(h,tpG,'edgecolor','b','linewidth',1.5)
% % highlight(h,fpG,'edgecolor','c')  % SUBGRAPH MUST HAVE THE SAME NODES AS THE GRAPH TO HIGHLIGHT
% % highlight(h,fnG,'edgecolor','m')
% % plot(fnG,'layout','force','iteration',365);
% highlight(h,fnG.Edges.EndNodes(:,1),fnG.Edges.EndNodes(:,2),'edgecolor','r','linestyle','--','linewidth',1)
% highlight(h,fpG.Edges.EndNodes(:,1),fpG.Edges.EndNodes(:,2),'edgecolor','g','linestyle','--','linewidth',1)

D=degree(GG);
[DD,IX]=sort(D);
highlight(h,IX(end-3:end),'nodecolor','m','markersize',6)


