function covar=covar(X,Y)
length_x=length(X);
X_bar=sum(X)/length_x;
Y_bar=sum(Y)/length_x;
covar=0;
for i=1:length(X)
    covar=covar+(1/(length_x-1))*(X(i)-X_bar)*(Y(i)-Y_bar);
end
end
