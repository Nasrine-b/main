% function qui calcule la variance d'un groupe de vecteurs
% par rapport à un vecteur moyen prédéterminé

function var = variance(vect, vmoy)
% vect : liste des angles des vecteurs
% vmoy : angle du vecteur moyen précédemment calculé 
%        (par ex. en passant par les coord cart)
% tout en degrés

    % on passe tous les angles dans ]vmoy-pi; vmoy+pi]
    vmoy = mod( vmoy+pi, 2*pi ) -pi; % ça devrait faire 0 mais pas tout à fait en pratique
    alpha = vmoy-pi;
    vect = mod( vect-alpha, 2*pi ) + vmoy-pi;
    % formule de la variance
    var = 1/length(vect) * sum( (vect-vmoy).^2 );

end