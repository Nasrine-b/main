%% Générer des réalisations de loi normale + plot 
% m=1;
% map(1,:)=[rand rand rand];
% colormap(map);
% sigma_m=1;
% nb=700;
% R=normrnd(m,sqrt(sigma_m),nb,15);
% hist(R)
% hold on
%x=m-4*sqrt(sigma_m):0.01:m+4*sqrt(sigma_m);
%p=nb*normpdf(x,m,sqrt(sigma_m));
%plot(x,p,'k');
%box off
%hold off
%% Produit cartésien 
V=[2,4]
V2=[1,2]
[X,Y]=meshgrid(V,V2)
res=[X(:), Y(:)]
%% Conversion des coordonnées pol->cart cart->pol
[m(1),m(2)]=cart2pol(V(1),V(2))
[P(1),P(2)]=pol2cart(m(1),m(2)) 
%% Covariance : utiliser les codes précédents
M=ones(4)
U=triu(M,1)