#include <iostream>
#include <vector>
#include <string>
#include "potential.h"
int main(int argc,char* argv[]){
	iteraction temp;
	int n=std::stoi(argv[1]);//iteraction with nthe layer. And without iteraction with n+1 layer.
	//let r go from 0 to r_cut to find this the minumum with this kind of iteraction
	std::vector<double> amplifier(n+2,0);
	std::vector<int> layernum(n+2,0);
	//store the coefficient for different scaled the lattice constant.
	for(size_t i=1;i<n+2;i++){
		temp=dist_n(2,i);
		amplifier[i]=temp.length;
		layernum[i]=temp.num;
	//	std::cout<<"neighbor layer: "<<i<<" numbers of neighbor: "<<temp.num<<" distance: "<<temp.length<<std::endl;
	}
	//started to find the mimimum value!!! very interesting.
	double la_const_equili=0;
	double la_energy=10000000;
	double r_stop=7.5;
	double r_start=0.0001;
	double r_step=0.00001;
	double r_temp=0;
	double sum_e=0;
	int N=int(r_stop-r_start)/r_step;
	for(int i=0;i<N;i++){
		r_temp=r_start+i*r_step;
		sum_e=0.0;
		for(int j=1;j<=n;j++){
			sum_e=sum_e+phi(amplifier[j]*r_temp/2)*layernum[j];
		}
		if(sum_e<la_energy){
			la_energy=sum_e;
			la_const_equili=r_temp;
		}
	}
	if(amplifier[n+1]*la_const_equili/2 > r_stop && amplifier[n]*la_const_equili/2 < r_stop)
	{
		std::cout<<"n layer iteration:"<<n<<" "<<"the lattice parameter is"<<la_const_equili<<std::endl;
	}
	else{
		std::cout<<"no solution for this situation: layer equals "<<n<<std::endl;
}
}
