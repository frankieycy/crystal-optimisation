#ifndef TOOL
#define TOOL

void line(){
	printf("----------------------------------------------\n");
}

void fline(FILE *f){
	fprintf(f,"----------------------------------------------\n");
}

double mod0(double x, double d){
	/* modulus mapped to [-d/2,d/2] */
	while(x>+d/2) x-=d;
	while(x<-d/2) x+=d;
	return x;
}

double min0(double a, double b){
	/* min of two num */
	if(a<b) return a;
	return b;
}

double roundsf(double a, int n){
	/* round to n sig fig */
	return roundf(a*pow(10,n))/pow(10,n);
}

double uniform(double a, double b){
	/* uniform sampling in [a,b] */
	return a+(b-a)*rand()/RAND_MAX;
}

#endif
