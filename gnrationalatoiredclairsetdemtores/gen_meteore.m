% génération de vecteurs de météores
% avec norme et angle du vecteur moyen
% et leur ecart-types

function mat = gen_meteore(norme, signorme, theta, sigtheta, nb_vect)
  
  mat = zeros(nb_vect, 2);
  for i=1:nb_vect
    mat(i,:) = [normrnd(norme, signorme) normrnd(theta, sigtheta)];
  end
 
 end