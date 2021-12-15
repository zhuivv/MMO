pmatrix = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/psysbp_10000.csv');
pmatrix(:,1) = [];

tmatrix = readmatrix('~/Desktop/mediation/DataMediationAnalysis/stdReHo/tmat.csv');
tmatrix(:,1) = [];
reho = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/ReHo.csv');
reho(:,1) = [];
cbf = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/CBF.csv');
cbf(:,1) = [];

ind = readmatrix('~/Desktop/mediation/DataMediationAnalysis/stdReHo/cvr_ind.csv');
ind(:,1) = [];

tfdr = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/realdataanlys/tfdr.csv');tfdr(:,1) = [];tfdr = tfdr';
sfdr = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/realdataanlys/sfdr.csv');sfdr(:,1) = [];sfdr = sfdr';

t3step = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/realdataanlys/t3step.csv');t3step(:,1) = [];t3step = t3step';
s3step = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/realdataanlys/s3step.csv');s3step(:,1) = [];s3step = s3step';

names = readmatrix('~/Dropbox/Mac/Desktop/mediation/updated/regionname.csv');
names(:,1) = [];

%% cor plots
figure; imagesc(corr(reho)); colormap parula;caxis([0 1]);colorbar;
figure; imagesc(corr(cbf)); colormap parula;caxis([0 1]);colorbar;

%% figure simulation
pmat= readmatrix('~/Dropbox/Mac/Desktop/mediation/simulation/Sep27/pmatrix_eff0.3_size15_N100/pmat2.csv');
pmat(:,1) = [];
figure; imagesc(pmat(datasample(1:100,100,'Replace',false),datasample(1:100,100,'Replace',false))); colormap parula;caxis([0 5]);colorbar;xlabel('Outcomes','Fontsize',28);ylabel('Mediators','Fontsize',28);

c_vec = 0.1:0.1:10;
r=-log(0.01);
lambda_vec = 1:0.1:5;
[s_rev,t_rev, c_rev, lambda_rev]=greedy_lik_fun(pmat,c_vec,lambda_vec,r,'KL');
s_in = setdiff(1:size(pmat,1),s_rev);
t_in = setdiff(1:size(pmat,2),t_rev);
figure; imagesc(pmat([s_in setdiff(1:100,s_in)],[t_in setdiff(1:100,t_in)])); colormap parula; caxis([0 5]);colorbar; xlabel('Outcomes','Fontsize',28);ylabel('Mediators','Fontsize',28);

%% real data set
W = pmatrix;
figure; imagesc(W); colormap parula;caxis([0 5]);colorbar;xlabel('ReHo','Fontsize',28); ylabel('CBF','Fontsize',28);
c_vec = 0.1:0.1:10;
r=-log(0.01);
lambda_vec = 1:0.1:5;
[s_rev,t_rev, c_rev, lambda_rev]=greedy_lik_fun(W,c_vec,lambda_vec,r,'KL');
s_in = setdiff(1:size(W,1),s_rev);
t_in = setdiff(1:size(W,2),t_rev);
%[s_out, t_out, c_max]=greedy_bipar(W,c_vec,10);
%s_in = setdiff(1:size(W,1),s_out);
%t_in = setdiff(1:size(W,2),t_out);
%figure;imagesc(W(s_in,t_in));colormap jet;colorbar;
% plot after shulffling
figure; imagesc(W([[s_in 63 93] sort(setdiff(1:107,[s_in 63 93]))],[t_in sort(t_rev)])); colormap parula; caxis([0 5]);xlabel('ReHo','Fontsize',28);ylabel('CBF','Fontsize',28); colorbar;  


figure; imagesc(W([sfdr sort(setdiff(1:107,sfdr))],[tfdr sort(setdiff(1:107,tfdr))])); colormap jet; caxis([0 9]);xlabel('ReHo','Fontsize',28);ylabel('CBF','Fontsize',28); colorbar;  

figure; imagesc(W([s3step sort(setdiff(1:107,s3step))],[t3step sort(setdiff(1:107,t3step))])); colormap jet; caxis([0 9]);xlabel 'ReHo';ylabel 'CBF'; colorbar;  


WW = W(s_in,t_in);
r=-log(0.02);
lambda_vec = 1:0.1:5;
[s_lasw,t_lasw, c_revw, lambda_lasw]=greedy_lik_fun(WW,c_vec,lambda_vec,r,'KL');
s_inld = setdiff(1:size(WW,1),s_lasw);
t_inld = setdiff(1:size(WW,2),t_lasw);
%figure;imagesc(WW(s_inw,t_inw));colormap jet;colorbar;
% plot after shulffling
figure; imagesc(WW([s_inld sort(s_lasw)],[t_inld sort(t_lasw)])); colormap jet; colorbar;


srev_ldfinal = setdiff(1:size(W,1),sin_ldfinal);
trev_ldfinal = setdiff(1:size(W,2),tin_ldfinal);
figure; imagesc(W([sin_ldfinal sort(srev_ldfinal)],[tin_ldfinal' sort(trev_ldfinal)])); colormap jet; caxis([0 4.5]); colorbar; 


writematrix(s_in,'~/Dropbox/Mac/Desktop/mediation/updated/sysbps.csv');
writematrix(t_in,'~/Dropbox/Mac/Desktop/mediation/updated/sysbpt.csv');


%% permutation
folder='~/Dropbox/Mac/Desktop/mediation/updated/seperate_anly/pmatf/';
folder2 = '~/Dropbox/Mac/Desktop/mediation/updated/seperate_anly/index/';
a=dir([folder 'pmat*.csv']);
total=size(a,1);

for i=1:total
pmatrix0 = readmatrix(fullfile(folder,sprintf('%s%i.csv','pmat',i)));
pmatrix0(:,1) = [];
%figure; imagesc(pmatrix0);

% select region
W = pmatrix0;
c_vec = 0.1:0.1:10;
r=-log(0.05);
lambda_vec = 1:0.1:4.5;
 try
  [s_rev,t_rev, c_rev, lambda_rev]=greedy_lik_fun(W,c_vec,lambda_vec,r,'KL');
  s_in = setdiff(1:size(W,1),s_rev);
  t_in = setdiff(1:size(W,2),t_rev);
%figure;imagesc(W(s_in,t_in));colormap jet;colorbar;
  writematrix(s_in, fullfile(folder2,sprintf('%s%i.csv','s_in',i)));
  writematrix(t_in, fullfile(folder2,sprintf('%s%i.csv','t_in',i)));
 catch ERM
    fprintf('%i%s\n',i,ERM.message);
 end
end



% visualize selected regions




