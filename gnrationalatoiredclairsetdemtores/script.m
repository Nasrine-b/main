%% Donnees 
ro=1;
%theta=[0,10,-10,80,90,100,170,180,190,260,270,280,350,360,370];
%theta=[10,13,11,12,11,10];
%theta=[0,10,-10,80];
theta=[45,135]
%% Conversion en rad + produit cartésien
[RO,THETA]=meshgrid(ro,theta);
THETA=THETA*pi/180;
V=[THETA(:),RO(:)];
%% Conversion en cartesien
X=1.*cos(THETA);
Y=1.*sin(THETA);
%% Calcul de V_moy
X_moy=mean(X);
Y_moy=mean(Y);
V_moy=[X_moy,Y_moy];
%% Conversion de V_moy en polaire
ro_moy=sqrt(X_moy.^2+Y_moy.^2);
theta_moy=2*atan(Y_moy./(X_moy+sqrt(X_moy.^2+Y_moy.^2)))
%[theta_moy,ro_moy]=cart2pol(X_moy,Y_moy);
%% Ecart type
theta_moy_pi=mod(theta_moy+pi,2*pi)-pi
theta_2pi=mod(THETA-theta_moy_pi+pi,2*pi)+theta_moy_pi-pi;

ecart_type=0;
for i=1:length(theta_2pi)
    diff=theta_2pi(i)-theta_moy_pi;
    ecart_type=ecart_type+(diff)^2;
end
ecart_type=sqrt(ecart_type/(length(theta_2pi)))
%% Plot des vecteurs
polar(THETA,RO,'o')
