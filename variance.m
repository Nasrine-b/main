% function qui calcule la variance d'un groupe de vecteurs
% par rapport à un vecteur moyen prédéterminé

function var = variance(vect, vmoy)
% vect : liste des angles des vecteurs
% vmoy : angle du vecteur moyen précédemment calculé 
%        (par ex. en passant par les coord cart)
% tout en degrés

    % on passe tous les angles dans ]vmoy-pi; vmoy+pi]
    alpha = vmoy-pi;
    vect = mod( vect-alpha, 2*pi ) + alpha;
    vmoy = mod( vmoy-alpha, 2*pi ) + alpha; % ça devrait faire 0 mais pas tout à fait en pratique
    
    % formule de la variance
    var = 1/length(vect) * sum( (vect-vmoy).^2 );

end