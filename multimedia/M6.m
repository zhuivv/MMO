function [s_in t_in eff pmat]=M6(x,M,Y)
% input x is the univariate exposure
% input M is the multivariate mediators
% input Y is the multivariate outcomes

% output M_in is the selected mediators
% output Y_in is the selected outcomes
% output eff is the estimated mediation effect
% output pmat is the mediation p-value matrix


% step 1
pmat = [];
for i=1:size(M,2)
    for j=1:size(Y,2)
        [t p] = sobel(x,M(:,i),Y(:,j));
        pmat(i,j) = p(2);
    end
end
s_in = [];t_in = [];
W = -log(pmat);
c_vec = 0.1:0.1:10;
r=0.1*(-log(0.02))+0.2*(-log(0.05))+0.5*(-log(0.1));
try
lambda_vec = 1:0.1:5;
[s_rev,t_rev, c_rev, lambda_rev]=greedy_lik_fun(W,c_vec,lambda_vec,r,'KL');
s_in = setdiff(1:size(W,1),s_rev);
t_in = setdiff(1:size(W,2),t_rev);
catch ERM
    fprintf(ERM.message);
end
if isempty(s_in)&&isempty(t_in)
    try
    [s_rev,t_rev, c_rev]=greedy_bipar(W,c_vec,1.2);
    s_in = setdiff(1:size(W,1),s_rev);
    t_in = setdiff(1:size(W,2),t_rev);
    catch ERM
    fprintf('%s\n',ERM.message);
    end
end
if isempty(s_in)&&isempty(t_in)
    msg = 'Error. no pattern extracted';
    error(msg)
end

% step 2
[coeffr,scorer] = pca(M(:,s_in));
Yr = Y(:,t_in);
effM = [];
for jj = 1:size(Yr,2)
aor = corr(x,scorer(:,1));
sea = sqrt((1-aor^2)/(length(x)-2));
mat01 = [ones(length(x),1) x];
mat02 = [ones(length(x),1) x];
mres = scorer(:,1)- mat01*(mat01'*mat01)^-1*mat01'*scorer(:,1);
yres = Yr(:,jj)-mat02*(mat02'*mat02)^-1*mat02'*Yr(:,jj);
bor = corr(mres,yres);
seb = sqrt((1-bor^2)/(length(mres)-2));
%bor = corr(scorer(:,1),Yr(:,jj));
effM(jj) = aor*bor;
tmp1 = bor^2*sea^2+aor^2*seb^2;
tmp2 = sea^2*seb^2;
se = sqrt(tmp1+tmp2);
end
effM6 = mean(effM);
eff = effM6;


