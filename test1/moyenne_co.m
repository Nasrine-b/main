function moy=moyenne_co(covar)
n=size(covar,1); %la matrice est carrée
moy=0;
for i=1:n
    for j=1:n
        moy=moy+covar(i,j);
    end
end
moy=moy/numel(covar);
end