% pour les tests

% % orage
norme = 10 ; sigma = 1;
nb_vect = 500;
orage = gen_cercle(norme, sigma, nb_vect);
%polar(orage)
[X,Y]= pol2cart(orage(:,2), orage(:,1));
orage_cartesien=[X Y];
subplot(1,2,1)
scatter(orage_cartesien(:,1), orage_cartesien(:,2))
title('Orage cartésien')
% météore
norme = 5; signorme = 1;
theta = 1.3 ; sigtheta = 0.3;
nb_vect = 1000;
meteore = gen_meteore(norme, signorme, theta, sigtheta, nb_vect);
[X,Y]= pol2cart(meteore(:,2), meteore(:,1))
meteore_cartesien=[X Y]
subplot(1,2,2)
scatter(meteore_cartesien(:,1), meteore_cartesien(:,2))
title('Meteor cartésien')
axis([0 5 0 15])

% fct qui calcule la moyenne et l'écart-type pour retrouver les bonnes moyennes
%[rho_moy,theta_moy]=vecteur_moyen(meteore(:,1),meteore(:,2))
meteor_cart_moy=[mean(X) mean(Y)]
[theta_moy,rho_moy]=cart2pol(meteor_cart_moy(1),meteor_cart_moy(2))
ecart_type_moy=sqrt(variance(rho_moy,theta_moy))