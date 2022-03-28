function [tvalue pvalue] = sobel(iv,mv,dv)
%medmdl = fitlm(iv,mv);
mat1 = [ones(length(iv),1) iv];
coef1 = (mat1'*mat1)^-1*mat1'*mv;
se1 = sqrt((mat1'*mat1)^-1*mat1'*var(mv)*((mat1'*mat1)^-1*mat1')');
%result = table2array(medmdl.Coefficients);
a = coef1(2);
sa = se1(2,2);
%fitmdl = fitlm([iv mv],dv);
%result2 = table2array(fitmdl.Coefficients);
mat2 = [ones(length(iv),1) iv mv];
coef2 = (mat2'*mat2)^-1*mat2'*dv;
se2 = sqrt((mat2'*mat2)^-1*mat2'*var(dv)*((mat2'*mat2)^-1*mat2')');
b = coef2(3);
sb = se2(3,3);
tmp1 = b^2*sa^2+a^2*sb^2;
tmp2 = sa^2*sb^2;
zsob = (a*b)/sqrt(tmp1);
zaro = (a*b)/sqrt(tmp1+tmp2);
psob = normcdf(-abs(zsob))*2;
paro = normcdf(-abs(zaro))*2;
tvalue = [zsob;zaro];
pvalue = [psob;paro];








