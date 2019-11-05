function [res1,res2]=construct_Vpa(nb_iter,sigma_m,sigma_theta,m,theta)
% Construction de m
length_sm=length(sigma_m);
M_m=zeros(length_sm,nb_iter,length(m));
for k=1:length(m)
    for i=1:length_sm
        while true
             l=normrnd(m(k),sigma_m(i),nb_iter,1)';
             if(l>0)
                 break;
             end
        end
           M_m(i,:,k)=l;
    end
end

% Construction de theta
length_st=length(sigma_theta);
M_t=zeros(length_st,nb_iter,length(theta));
for k=1:length(theta)
    for i=1:length_st
       M_t(i,:,k)=normrnd(theta(k),sigma_theta(i),nb_iter,1)';
    end
end
res1=M_m;
res2=M_t;
end