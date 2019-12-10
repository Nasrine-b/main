% génération de vecteurs aléatoires en cercle 
% à partir d'une norme et un écart-type de norme 

function mat = gen_cercle(norme, sigma, nb_vect)
  
  mat = zeros(nb_vect, 2);
  for i=1:nb_vect
    mat(i,:) = [normrnd(norme, sigma), unifrnd(0, 2*pi())];
  end 
  
end 