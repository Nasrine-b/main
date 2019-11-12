ro=1;
theta=[0,10,-10,80,90,100,170,180,190,260,270,280,350,360,370];
[RO,THETA]=meshgrid(ro,theta);
THETA=THETA*pi/180;
V=[THETA(:),RO(:)];
[X,Y]=pol2cart(V(:,2),V(:,1));
X_moy=mean(X);
Y_moy=mean(Y);
V_moy=[X_moy,Y_moy];
[theta_moy,ro_moy]=cart2pol(X_moy,Y_moy);

theta_2pi=mod(THETA,2*pi)-pi;
theta_moy_pi=mod(theta_moy,2*pi)-pi;

diff=zeros(length(theta),1);
ecart_type=0;
for i=1:length(theta_2pi)
    diff(i)=theta_2pi(i)-theta_moy_pi;
    ecart_type=ecart_type+(diff(i))^2;
end
ecart_type=ecart_type/(length(theta_2pi)-1)

