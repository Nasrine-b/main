#include "parametre.h"

float f_S(float* vec,int taille,int p){
	float Sp=0;
	for (int i = 0; i < taille; ++i)
	{
		float temp=p*vec[i];
		Sp=Sp+sin(temp);
	}
	return Sp;
}

float f_C(float* vec,int taille,int p){
	float Cp=0;
	for (int i = 0; i < taille; ++i)
	{
		Cp=Cp+cos(p*vec[i]);
	}
	return Cp;
}

float f_T(float Cp_b, float Sp_b){
	float Tp=atan(Sp_b/Cp_b);

	if (Cp_b<0)
	{
		Tp=Tp+3.14;
	}
	else if(Cp_b>0 && Sp_b<0)
	{
		Tp=Tp+2*3.14;
	}

	return Tp;
}

float f_R(float C, float S){
	return sqrt(pow(C,2)+pow(S,2));
};

float f_V(float R){
	return 1-R;
}

float f_sigma(float R){
	return sqrt(-2*logf(R));
}

float f_delta(float T, float R){
	return (1-T)/(2*R*R);
}

float f_s(float R1, float R2, float T1, float T2){
	return (R2*sin(T2-2*T1))/pow((1-R1),3/2);
}

float f_k(float R1, float R2, float T1, float T2){
	return (R2*cos(T2-2*T1)-pow(R1,4))/pow((1-R1),3/2);
}

void circ_data(float* vec, int taille, float* R1, float* R1_b, float* theta, float* V, float* sigma, float* delta, float* s, float* k){

	//Creation Cp, Cp_b, Sp, Sp_b 
	float C1=f_C(vec,taille,1);
	float C2=f_C(vec,taille,2);
	float C1_b=C1/(float)taille;
	float C2_b=C2/(float)taille;

	float S1=f_S(vec,taille,1);
	float S2=f_S(vec,taille,2);
	float S1_b=S1/(float)taille;
	float S2_b=S2/(float)taille;

	float T1=f_T(C1_b,S1_b);
	float T2=f_T(C2_b,S2_b);

	//Calculs des retours
	*(R1)=f_R(C1,S1);
	float R2=f_R(C2,S2);

	*(R1_b)=*(R1)/(float)taille;
	float R2_b=R2/(float)taille;	

	*(theta)=T1;

	*(V)=f_V(*(R1_b));

	*(sigma)=f_sigma(*(R1_b));

	*(delta)=f_delta(T2,*(R1_b));

	*(s)=f_s(*(R1_b),R2_b,T1,T2);

	*(k)=f_k(*(R1_b),R2_b,T1,T2);
};