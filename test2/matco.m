function matco=matco(X,Y)
    matco=[covar(X,X),covar(X,Y);covar(Y,X),covar(Y,Y)];
end