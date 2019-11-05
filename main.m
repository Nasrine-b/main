% Donn�es 
m=1;
sigma_m=[1,2,3];
theta=[0,10,-10,80,90,100,170,180,190,260,270,280,350,360,370];
sigma_theta=[1,2,5,10,45,90,180,360];
nb_valeurs=[10,50,100,200];
%% Construction des Vp 
[X,Y]=meshgrid(m,theta);
Vp=[X(:),Y(:)]; 
%% Construction des Vpa pour un nb_valeurs
%length(sigma_m)
V=VectPol(m,10,360,theta(3),theta(length(theta)));
%[X,Y]=meshgrid(x,y);
%Vpa_pol=[X(:),Y(:)];
figure(1)
polar(V(2,:),V(1,:),'o')
%% Conversion Vpa_pol � Vpa_r
% [x,y]=pol2cart(Vpa_pol(:,1),Vpa_pol(:,2));
% Vpa_r=[x,y];
% figure(2)
% scatter(Vpa_r(:,1),Vpa_r(:,2))
% title('Dispersion des vercteurs vitesses')
%  %% Calcul des matrices de covariance 
% % M=matrice(Vpa_r);
% % M_moy=matrice_moyenne_co(M)
% %% Calcul d'une matrice de covariance 
% X=Vpa_r(:,1);
% Y=Vpa_r(:,2);
% M=matco(X,Y);
% Moy_M=Mean_matrix(M); %ne sert a rien
% %% Coefficient de correlation 
% coef_cor=M(1,2)/(sqrt(M(1,1))*sqrt(M(2,2))) %proche de 0 : correlation faible 

