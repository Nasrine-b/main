#include "parametre.h"

int main(int argc, char const *argv[])
{
	int taille=5;
	float* vec=malloc(taille*sizeof(int));
	vec[0]=40;
	vec[1]=20;
	vec[2]=35;
	vec[3]=60;
	vec[4]=12;

	float R1, R1_b, theta, V, sigma, delta, s, k;
	float* R1_pt, * R1_b_pt, * theta_pt, * V_pt, * sigma_pt, * delta_pt, * s_pt, * k_pt;

	R1_pt=&(R1);
	R1_b_pt=&(R1_b);
	theta_pt=&(theta);
	V_pt=&(V);
	sigma_pt=&(sigma);
	delta_pt=&(delta);
	s_pt=&(s);
	k_pt=&(k);
	
	circ_data(vec,taille,R1_pt,R1_b_pt,theta_pt,V_pt,sigma_pt,delta_pt,s_pt,k_pt);

	printf("R1=%f, R1_b=%f, theta=%f, V=%f, sigma=%f, delta=%f, s=%f, k=%f \n",R1,R1_b,theta,V,sigma,delta,s,k);

	return 0;
}