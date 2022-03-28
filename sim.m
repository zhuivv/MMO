function [x M Y] = sim(n,a,b,c,c1,c2,cat)
% input: n is the number of subjects
% input: a is the partial correlation between M and x
% input: b is the partial correlation between M and Y
% input: c is the partial correlation between x and Y
% input: c1,c2 are the coefficients for first and third block
% input: cat is the cardinality size

% output: x is the simulated univariate exposure
% output: M is the simulated multivariate mediators
% output: Y is the simulated multivariate outcomes

Omg=eye(3)+squareform([a,c,b]);
D = mvnrnd([0 0 0 ],Omg^-1,n);
m1 = mvnrnd(0.1,1,n);
m2 = mvnrnd(-0.1,1,n);
m0 = D(:,2);
y1 = mvnrnd(0,1,n);
y2 = mvnrnd(0,1,n);
y0 = D(:,3);

coef = (mvnrnd([0 0 0],eye(3)*0.01,100))';
coef(1,1:30) = coef(1,1:30)+1; 
coef(2,31:(30+cat)) = coef(2,31:(30+cat))+1;
coef(3,(30+cat+1):100) = coef(3,(30+cat+1):100)+0.9;
M = [m1 m0 m2]*coef;


Y = [y1 y0 y2]*coef;
x = D(:,1);