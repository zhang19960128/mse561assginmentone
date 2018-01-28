#ifndef potential_h
#define potential_h
typedef struct iteraction{
	double length;
	int num;
}iteraction;//this is the datastructure to store the iteracting atoms and numbers of iteraction.
double phi(double);
iteraction dist_n(double,int);
#endif
