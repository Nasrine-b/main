% fichier de tests

% data
rho=1;
theta=[0,10,-10,80,90,100,170,180,190,260,270,280,350,360,390];
% conversion en radians de theta
theta = theta*pi/180;

% produit cartesien de toutes les associations (rho, theta)
[rho, theta] = meshgrid(rho, theta);
% affichage des vecteurs
polar(theta, rho, 'o')


% calcul du vecteur moyen
[rho_moy, theta_moy] = vecteur_moyen(rho, theta);

% calcul de l'écart-type des angles des vecteurs
sd = sqrt( variance(theta, theta_moy) )