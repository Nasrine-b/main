function [res]=matrice_moyenne_co(matrice)
res=0;
N=length(matrice);
res=zeros(N);
for i=1:N
    for j=1:N
           % res=res+matrice(1,2,i,j)/N^2;
           res(i,j)=Mean_matrix(matrice(:,:,i,j));
    end
end
end