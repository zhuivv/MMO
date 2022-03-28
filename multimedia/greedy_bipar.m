function [s_out t_out c_max]=greedy_bipar(W,c_vec,lambda)
W1 = W;
clear remove_s_vec remove_t_vec dST_vec
for l=1:size(c_vec,2)
    l;
    W=abs(W1);
    c=c_vec(l);
    density_list = zeros((size(W,1)+size(W,2))-1,8);
    for i=1:(size(W,1)+size(W,2))
        if W==0
            break
        end
        % mean of columns
        C=sum(W,1);
        % mean of rows
        R=sum(W,2);
        %[dT,IndT]=min(C(C>0));
        % min sum in column 
        dT=min(C(C>0));
        IndT=find(C==dT);
        % min sum in row 
        dS=min(R(R>0));
        IndS=find(R==dS);
        %c=|S|/|T|=30/50
        %c=1;%sum(R~=0)/sum(C~=0);
        if sqrt(c)*dS <= dT/sqrt(c)
            W(IndS(1),:)=zeros(1,size(W,2));
            density_list(i,1)=sum(sum(W,2)~=0);
            density_list(i,2)=sum(sum(W,1)~=0);
            density_list(i,3)=1;
        else
            W(:,IndT(1))=zeros(size(W,1),1);
            density_list(i,1)=sum(sum(W,2)~=0);
            density_list(i,2)=sum(sum(W,1)~=0);
            density_list(i,3)=2;
        end
        density_list(i,4)=IndS(1);
        density_list(i,5)=dS;
        density_list(i,6)=IndT(1);
        density_list(i,7)=dT;
        density_list(i,8)=(sum(sum(W,2)))/(sqrt(density_list(i,1)*density_list(i,2)))^lambda;
    end


    %ouput W
    [dST, indST]=max(density_list(:,8));
    remove_s=[];
    remove_t=[];
    for j=1:indST
        if density_list(j,3)==1
            remove_s =[remove_s density_list(j,4)];
        else
            remove_t =[remove_t density_list(j,6)];
        end
end
remove_s_vec{l}=remove_s;
remove_t_vec{l}=remove_t;
dST_vec(l)=dST;
end
c_max = find(dST_vec == max(dST_vec));
s_out = remove_s_vec{c_max(1)};
t_out = remove_t_vec{c_max(1)};
end
