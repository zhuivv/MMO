function [p permt t0] =permutation(x,M,Y,B);
% B is the permutation repetations
% output p is the permutation p value
% permt is the permutation t statistic
% t0 is the t statistic from the original data


% original test statistic
[M_in Y_in eff pmat] = M6(x,M,Y);
t0 = abs(sum(-log(pmat(M_in,Y_in)),[1 2])/sqrt(length(M_in)*length(Y_in)));


% permutation
permt=[];
parfor i=1:B
Mind = datasample(1:size(M,1),size(M,1));
xind = datasample(1:length(x),length(x));

M_perm = M(Mind,:); x_perm = x(xind);
[M_pin Y_pin effp pmatp] = M6(x_perm,M_perm,Y);
permt(i) = abs(sum(-log(pmatp(M_pin,Y_pin)),[1 2])/sqrt(length(M_pin)*length(Y_pin)));

end

p = sum(permt>t0)/B;




























