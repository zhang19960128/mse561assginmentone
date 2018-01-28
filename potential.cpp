#include <iostream>
#include <vector>
#include "potential.h"
#include <cmath>
#include <algorithm>
double phi(double r){
	double r0=7.0;
	double r_cut=7.5;
	double sigma=3.405;
	double epsilon=0.010323;
	double A=-6.8102128*1e-3;
	double B=-5.5640876*1e-3;
	if(r<r0){
		return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
	}
	else if(r>r_cut){
		return 0;
	}
	else{
		return A*pow(r-r_cut,3)+B*pow(r-r_cut,2);
	}
}
iteraction dist_n(double rinside,int n){
//determine the lenght of the n layer iteracting atoms and numbers of interacting. rinside is the crystal constant
	double r=rinside/2;
	//the whole configuration has two different slice.
	double onelayer[2*n][2*n][2];
	double onelen[2*n][2*n];
	double twolayer[2*n][2*n][2];
	double twolen[2*n][2*n];
	for(int i=-n;i<n;i++)
		for(int j=-n;j<n;j++){
		onelayer[i+n][j+n][0]=2*i*r;
		onelayer[i+n][j+n][1]=2*sqrt(3)*j*r;
		twolayer[i+n][j+n][0]=(2*i+1)*r;
		twolayer[i+n][j+n][1]=((2*j+1)*sqrt(3))*r;
		}
	for(int i=-n;i<n;i++)
		for(int j=-n;j<n;j++){
			onelen[i+n][j+n]=onelayer[i+n][j+n][0]*onelayer[i+n][j+n][0]+onelayer[i+n][j+n][1]*onelayer[i+n][j+n][1];
			twolen[i+n][j+n]=twolayer[i+n][j+n][0]*twolayer[i+n][j+n][0]+twolayer[i+n][j+n][1]*twolayer[i+n][j+n][1];
		}
	//we have collected all the information about different distance according to different lenght;
	std::vector<double> all(2*2*n*2*n,0);
	for(size_t i=0;i<2*n;i++)
		for(size_t j=0;j<2*n;j++){
			all[i*2*n+j]=onelen[i][j];
			all[i*2*n+j+2*n*2*n]=twolen[i][j];
		}
	std::sort(all.begin(),all.end());
	double len=1;
	int count=0;
	int layer=0;
	int record_num=0;
	int record_len=0;
	for(std::vector<double>::iterator a=all.begin()+1;a!=all.end();a++){
		if(fabs(*a-len)>1e-10){
			layer++;
			count=1;
			len=*a;
		}
		else{
			count++;
			len=*a;
		}
		if(layer==n+1) break;
		record_num=count;
		record_len=len;
	}
	if(n==1){
		record_len=(2*r)*(2*r);
		record_num=6;
	}
	iteraction temp;
	temp.length=sqrt(record_len);
	temp.num=record_num;
	return temp;
}
