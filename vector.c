#ifndef VECTOR
#define VECTOR

/*
Note:
	- vec: dynamic double array
*/

typedef struct vec{
	int l;
	double *v;
} vec;

vec v(int l){
	/* init empty vec */
	vec x={l,malloc(sizeof(double)*l)};
	return x;
}

vec vv(int l, double *v){
	/* init filled vec */
	vec x={l,v};
	return x;
}

/*************/

/* arrays: 2d,vec,etc. */

int **int2d(int r, int c){ 
	/* init 2d int array */
	int **x;
	x=malloc(sizeof(int*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(int)*c);
	return x;
}

double *double1d(int r){
	/* init 1d double array */
	double *x;
	x=malloc(sizeof(double)*r);
	return x;
}

double **double2d(int r, int c){ 
	/* init 2d double array */
	double **x;
	x=malloc(sizeof(double*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(double)*c);
	return x;
}

vec *vec1d(int l){
	/* init 1d vec array */
	vec *x;
	x=malloc(sizeof(vec)*l);
	return x;
}

vec **vec2d(int r, int c){
	/* init 2d vec array */
	vec **x;
	x=malloc(sizeof(vec*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(vec)*c);
	return x;
}

/*************/

/* handling */

void freev(vec a){free(a.v);}

vec cart(double a, double b, double c){ 
	/* Cartesian vec */
	vec x=v(3);
	x.v[0]=a; x.v[1]=b; x.v[2]=c;
	return x;
}

vec z(int l){
	/* zero vec */
	vec x=v(l);
	for(int i=0; i<l; i++) x.v[i]=0;
	return x;
}

vec randv(int l, double eps){
	/* random vec */
	vec x=v(l);
	for(int i=0; i<l; i++) x.v[i]=uniform(-eps,eps);
	return x;
}

vec copy(vec a){
	/* copied vec */
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i];
	return x;
}

/*************/

/* arithmetic */

double min(vec a){
	/* min entry */
	double m=a.v[0];
	for(int i=0; i<a.l; i++)
		if(a.v[i]<m) m=a.v[i];
	return m;
}

double mag(vec a){
	/* vec norm */
	double l2=0;
	for(int i=0; i<a.l; i++) l2 += a.v[i]*a.v[i];
	return sqrt(l2);
}

vec add(vec a, vec b){
	/* add vec */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in add(vec a, vec b): a and b not of same dimension.\n");
	}
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]+b.v[i];
	return x;
}

vec minus(vec a, vec b){
	/* subtract vec */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in minus(vec a, vec b): a and b not of same dimension.\n");
	}
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]-b.v[i];
	return x;
}

vec vmul(vec a, double b){
	/* scalar mult */
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]*b;
	return x;
}

vec vdiv(vec a, double b){
	/* scalar div */
	return vmul(a,1/b);
}

vec unit(vec a){
	/* unit vec */
	double l=mag(a);
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]/l;
	return x;
}

double dot(vec a, vec b){
	/* dot prod */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in dot(vec a, vec b): a and b dimensions not matched.\n");
		exit(-1);
	}
	double p=0;
	for(int i=0; i<a.l; i++) p += a.v[i]*b.v[i];
	return p;
}

vec cross(vec a, vec b){
	/* cross prod */
	if(a.l!=3||b.l!=3){
		fprintf(stderr,
		"Error in cross(vec a, vec b): a or b not of dimension 3.\n");
		exit(-1);
	}
	vec x=v(a.l);
	x.v[0]=a.v[1]*b.v[2]-a.v[2]*b.v[1];
	x.v[1]=a.v[2]*b.v[0]-a.v[0]*b.v[2];
	x.v[2]=a.v[0]*b.v[1]-a.v[1]*b.v[0];
	return x;
}

double triple(vec a, vec b, vec c){
	/* triple prod */
	if(a.l!=3||b.l!=3||c.l!=3){
		fprintf(stderr,
		"Error in triple(vec a, vec b, vec c): a, b or c not of dimension 3.\n");
		exit(-1);
	}
	double v=dot(cross(a,b),c);
	return v;
}

double vcos(vec a, vec b){
	/* cos angle between */
	return dot(a,b)/(mag(a)*mag(b));
}

void printv(vec a){
	printf("{");
	for(int i=0; i<a.l; i++){
		printf(" %lf ",a.v[i]);
		if(i<a.l-1) printf(",");
	}
	printf("}\n");
}

#endif
