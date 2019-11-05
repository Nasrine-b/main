m=6;
sigma=1;
theta=[0,10,-10,80,90,100,170,180,190,260,270,280,350,360,370];

VectsPol=zeros(2,360);

for i=1:360
    VectsPol(1,i)=normrnd(m,sigma);
    VectsPol(2,i)=unifrnd(-10,370);
end

polar(VectsPol(2,:),VectsPol(1,:),'o')