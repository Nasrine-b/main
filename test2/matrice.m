function matrice=matrice(vecteurs)
matrice=zeros(2,2,length(vecteurs),length(vecteurs));
n=length(vecteurs);
for i=1:n-1
   for j=i+1:n
       vect1=vecteurs(i,:);
       vect2=vecteurs(j,:);
       mat=matco(vect1,vect2);
       matrice(:,:,i,j)=mat;
       matrice(:,:,j,i)=mat;
   end
end
for i=1:n
       vect1=vecteurs(i,:);
       mat=matco(vect1,vect1);
       matrice(:,:,i,i)=mat;
end
end