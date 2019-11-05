function [res]=VectPol(m,sigma_m,taille,theta_min,theta_max)
res=zeros(2,taille);

for i=1:taille
    res(1,i)=normrnd(m,sigma_m);
    res(2,i)=unifrnd(theta_min,theta_max);
end
end