function [s_rev,t_rev, c_rev, lambda_rev]=greedy_lik_fun(W,c_vec,lambda_vec,r,method)
% W is the input matrix can be adjacency matrix or correlation matrix;
% c_vec is the values for c in greedy search;
% lambda_vec is the lambda to tune;
% r is the threshold;
% method can take 'LogLik', 'entropy', 'KL';
W1 = W;
N=size(W,1);
M=size(W,2);
E=ones(size(W,1),size(W,2));
clear remove_s_cvec remove_t_cvec cmax_vec lik_vec
for k=1:size(lambda_vec,2)
    W=abs(W1);
    lambda=lambda_vec(k);
    %c_max
    [s_out, t_out, c_max]=greedy_bipar(W,c_vec,lambda);
    %likelihood
    s_in = setdiff(1:N,s_out);
    t_in = setdiff(1:M,t_out);
    %A
    A=W>r;
    At=1-A;
    %pi_MLE
    pi_s=sum(sum(A(s_in,t_in),2))/sum(sum(E(s_in,t_in),2));
    pi_0=(sum(sum(A,2))-sum(sum(A(s_in,t_in),2)))/(sum(sum(E,2))-sum(sum(E(s_in,t_in),2)));
    p = sum(sum(A,2))/sum(sum(E,2));
    %L = pi_s^(sum(sum(A(s_in,t_in),2)))*(1-pi_s)^(sum(sum(At(s_in,t_in),2)))*pi_0^(sum(sum(A(s_out,t_out),2)))*(1-pi_0)^(sum(sum(At(s_out,t_out),2)));
    switch method
        case 'LogLik'
            logL = (sum(sum(A(s_in,t_in),2)))*log(pi_s)+(sum(sum(At(s_in,t_in),2)))*log(1-pi_s)+(sum(sum(A(s_out,t_out),2)))*log(pi_0)+(sum(sum(At(s_out,t_out),2)))*log(1-pi_0);
        case 'entropy'
            logL = length(s_in)*length(t_in)*(pi_s*log(pi_s)+(1-pi_s)*log(1-pi_s))+length(s_out)*length(t_out)*(pi_0*log(pi_0)+(1-pi_0)*log(1-pi_0));
        case 'KL'
            logL = (sum(sum(A(s_in,t_in),2)))*(pi_s*log(pi_s/p))+(sum(sum(At(s_in,t_in),2)))*(1-pi_s)*log((1-pi_s)/(1-p))+(sum(sum(A(s_out,t_out),2)))*(pi_0*log(pi_0/p))+(sum(sum(At(s_out,t_out),2)))*(1-pi_0)*log((1-pi_0)/(1-p));
        %case 'KLq'
         %   logL = -((sum(sum(A(s_in,t_in),2)))*(p*log(p/pi_s))+(sum(sum(At(s_in,t_in),2)))*(1-p)*log((1-p)/(1-pi_s))+(sum(sum(A(s_out,t_out),2)))*(p*log(p/pi_0))+(sum(sum(At(s_out,t_out),2)))*(1-p)*log((1-p)/(1-pi_0)));
        otherwise 
            error('method: Wrong input');
    end
    %output_a
    %if isnan(logL)
    %    logL = 100000; 
    %end
    remove_s_cvec{k}=s_out;
    remove_t_cvec{k}=t_out;
    cmax_vec{k}=c_vec(c_max);
    lik_vec(k)=logL;
end
    %lambda_list=repelem(lambda_vec,size(c_vec,2));
    %c_list=repmat(c_vec,1,size(lambda_vec,2));
    %output_b
    lambda_idx=find(lik_vec == max(lik_vec));
    s_rev= remove_s_cvec{lambda_idx(1)};
    t_rev = remove_t_cvec{lambda_idx(1)};
    c_rev = cmax_vec{lambda_idx(1)};
    lambda_rev = lambda_vec(lambda_idx);
end







