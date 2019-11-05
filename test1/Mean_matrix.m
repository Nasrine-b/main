function res=Mean_matrix(M)
N=length(M);
res=0;
for i=1:N
    for j=1:N
        res=res+M(i,j);
    end
end
res=res/N^2;
end