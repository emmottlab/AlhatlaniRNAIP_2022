% Bader Corona RNA IP

clear
clc

dat = struct();
path = '/Users/ed/Dropbox/Liverpool/Collaborations/Bader_RNAIP/txt/';

% Load data

dat.rawprot = readtable([path , 'proteinGroups.txt']);
dat.prot = dat.rawprot;
%% Data cleanup
% Data were searched with 16 TMT channels, only 9 were used:
% 3 each of Un (control), M (MERS) and S (SARS-CoV-2)

% QC
temp = table2array(dat.prot(:,19:34));
figure
bar(log10(sum(temp)));

clear temp
%%
% Out of 1:16: I want: [10,6,13,15,9,5,16,11,7]
TMTorder = [10,6,13,15,9,5,16,11,7];

% Remove Contaminants and Reverse
dat.prot.Reverse = categorical(dat.prot.Reverse);
dat.prot.PotentialContaminant = categorical(dat.prot.PotentialContaminant);

dat.prot = dat.prot(dat.prot.Reverse ~= '+',:);
dat.prot = dat.prot(dat.prot.PotentialContaminant ~= '+',:);

dat.mat = table2array(dat.prot(:,19:34));

% Now reorder to account for TMT randomisation
dat.mat = dat.mat(:,TMTorder);

%% Data normalisation
dat.mat(dat.mat == 0) = NaN;

% Column normalise by dividing by the median
dat.mat = dat.mat ./ nanmedian(dat.mat);

% Remove rows with 6+ NaN;
logNaN = sum(isnan(dat.mat),2) >= 6;

% Remove rows
dat.mat = dat.mat(~logNaN,:);
dat.prot = dat.prot(~logNaN,:);

% KNNimpute
dat.mat = knnimpute(dat.mat);

% Sample ID
mock = [1,2,4];
mers = [3,5,6];
sars = [7:9]; 

% Normalise by dividing by the control
dat.mat = dat.mat ./ mean(dat.mat(:,mock),2,'omitnan');

% Log2
dat.mat = log2(dat.mat);

%% QC
[coeff,score,latent,tsquared,explained,mu] = pca(dat.mat);


% figure
% scatter(coeff(mock,1),coeff(mock,2),'filled','k')
% hold on
% 
% scatter(coeff(mers,1),coeff(mers,2),'filled','b')
% 
% scatter(coeff(sars,1),coeff(sars,2),'filled','r')
% 
% xlabel('PCA 1')
% ylabel('PCA 2')

%% Stats
[h,p,ci,stats] = ttest2(dat.mat(:,mock), dat.mat(:,mers),'Dim',2);
[fdr,q] = mafdr(p);

figure


subplot(2,2,1)
scatter(mean(dat.mat(:,mers),2),-log10(q),'filled','k','MarkerFaceAlpha',0.3)
title('MERS over Mock')
xlabel('Log_2 MERS over Mock')
ylabel('-log_1_0 q-value')
xlim([-4,4])

[h2,p2,ci2,stats2] = ttest2(dat.mat(:,mock), dat.mat(:,sars),'Dim',2);
[fdr2,q2] = mafdr(p2);

subplot(2,2,2)
scatter(mean(dat.mat(:,mers),2),-log10(q2),'filled','k','MarkerFaceAlpha',0.3)
title('SARS-CoV-2 over Mock')
xlabel('Log_2 SARS-CoV-2 over Mock')
ylabel('-log_1_0 q-value')
xlim([-4,4])

subplot(2,2,[3:4])
scatter(coeff(mock,1),coeff(mock,2),'filled','k')
hold on

scatter(coeff(mers,1),coeff(mers,2),'filled','b')

scatter(coeff(sars,1),coeff(sars,2),'filled','r')
hold off
xlabel('PCA 1')
ylabel('PCA 2')
title('Quality control: PCA analysis')
legend({'Mock','MERS','SARS-CoV-2'})

%% Figs for Bader:
% Add gene name annotations
dat.prot.GN = extractAfter(dat.prot.FastaHeaders,'GN=');
dat.prot.GN = extractBefore(dat.prot.GN,' PE=');


% 1. Better volcano plot

[h,p,ci,stats] = ttest2(dat.mat(:,mock), dat.mat(:,mers),'Dim',2);
[fdr,q] = mafdr(p);

fcCutoff = 0.5;
qCutoff = 0.05;

lgFcCut = mean(dat.mat(:,mers),2) >= fcCutoff;
lgQCut = q < qCutoff;
cutoff = lgFcCut + lgQCut == 2;

fcCutoff2 = 1;
lgFcCut2 = mean(dat.mat(:,mers),2) >= fcCutoff2;
txtcutoff = lgFcCut2 + lgQCut == 2;

figure
scatter(mean(dat.mat(~cutoff,mers),2),-log10(q(~cutoff)),'filled','k','MarkerFaceAlpha',0.5)
hold on
scatter(mean(dat.mat(cutoff,mers),2),-log10(q(cutoff)),'filled','r','MarkerFaceAlpha',0.5)

line([fcCutoff,fcCutoff],[0,3],'LineStyle','--','Color','k')
line([-fcCutoff,-fcCutoff],[0,3],'LineStyle','--','Color','k')
line([-4,4],[-log10(qCutoff),-log10(qCutoff)],'LineStyle','--','Color','k')
text(mean(dat.mat(txtcutoff,mers),2),-log10(q(txtcutoff)),dat.prot.GN(txtcutoff));
title('MERS 5'' RNA over Control')
xlabel('Mean Log_2 MERS over Control')
ylabel('-log_1_0 q-value')
xlim([-3,3])
ylim([0 3])
set(gca,'FontSize',14)
% 2. Tables of high confidence hits

%% Export table
dat.prot.MeanMERSoverMock = mean(dat.mat(:,mers),2);
dat.prot.ttestQvalue = q;
writetable(dat.prot(cutoff,:),'/Users/ed/Dropbox/Liverpool/Collaborations/Bader_RNAIP/hits.csv')