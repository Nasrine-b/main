% fonction qui calcule le vecteur moyen d'un ensemble de vecteur polaires

function [rho_moy, theta_moy] = vecteur_moyen(rho, theta)
% rho et theta : tableaux des coordonées en polaire

    % passage en cartésien
    [x, y] = pol2cart(theta, rho);
        % formules : x = rho*cos(theta)
        %            y = rho*sin(theta)

    % calcul du vecteur moyen 
    xm = mean(x);
    ym = mean(y);

    % retour en polaire
    [theta_moy, rho_moy] = cart2pol(xm, ym);
        % formules : rho_moy = sqrt(xm.^2 + ym.^2);
        %            theta_moy = 2* atan( ym./(xm + sqrt(xm.^2 + ym.^2)) )
    
end