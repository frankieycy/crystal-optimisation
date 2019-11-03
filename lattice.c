#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tool.c"
#include "vector.c"

/*
# Note
steepest descent (SD) & conjugate gradient (CG)
Lennard-Jones (LJ) & Tersoff (TE) potential
optimise: (1) fcc Ar (LJ)
          (2) dia Si (TE)

# Definitions
A    atom name (Ar,Si) (for record only)
S    crystal system name (sc,bcc,fcc,dia)
N    total number of atoms
R    atom coordinates (vec[N])
a    lattice constant (unit:angstrom)

K    nbor cutoff distance (unit:angstrom)
nb   num of nbor within cutoff
Nb   nbor indices (int[N][N])

V    pair pbc displacements (vec[N][N])
D    pair pbc distances (double[N][N])

# Test systems lookup
(1) LJ : Ar fcc 555 ; K=8.5   ; a=5.27
(2) TE : Si dia 333 ; K=999.0 ; a=5.43

search " [##] " to look for changeable params
*/

// periodicity [##]
#define Nx 3
#define Ny 3
#define Nz 3
#define a  5.43
#define K  999.

char *S,*A="X";
int  N,nb,**Nb;

double
// lattice vectors
t1[]={a*Nx,0,0},
t2[]={0,a*Ny,0},
t3[]={0,0,a*Nz},
omega,**D;

vec
a1,a2,a3, // direct
b1,b2,b3, // reciprocal
*R,**V;

/***************************************/

/* print */

void print_latt(FILE *f, char *cm){
	/* print coord in xyz format */
	fprintf(f,"%d\n",N);
	fprintf(f,"%s\n",cm);
	for(int i=0; i<N; i++)
		fprintf(f,"H %lf %lf %lf\n",
		R[i].v[0],R[i].v[1],R[i].v[2]);
}

void print_nbor(char *cm){
	/* print neighbor list */
	FILE *f=fopen("nbor.txt","w");
	fprintf(f,"# %s %s %d%d%d; K=%.4f; %s\n",A,S,Nx,Ny,Nz,K,cm);
	for(int i=0; i<N; i++) // each atom
		for(int j=0; j<nb; j++) // its nbor
			fprintf(f,"%d %d %.4f\n",i,Nb[i][j],D[i][Nb[i][j]]);
	fclose(f);
}

void info(){
	/* print sys info */
	printf("A: %s\n",A);
	printf("S: %s\n",S);
	printf("a: %.2f\n",a);
	printf("N: %d%d%d [%d]\n",Nx,Ny,Nz,N);
}

/***************************************/

/* latt vectors */

void set_recip(){
	// omit 2π
	omega=triple(a1,a2,a3);
	b1=vdiv(cross(a2,a3),omega);
	b2=vdiv(cross(a3,a1),omega);
	b3=vdiv(cross(a1,a2),omega);
}

void init(){
	/* initialise direct space */
	a1=vv(3,t1);
	a2=vv(3,t2);
	a3=vv(3,t3);
	b1=v(3);
	b2=v(3);
	b3=v(3);
	set_recip();

	Nb=int2d(N,N);
	D=double2d(N,N);
	V=vec2d(N,N);
}

void print_lattv(){
	/* lattice vectors */
	printf("a1 = "); printv(a1);
	printf("a2 = "); printv(a2);
	printf("a3 = "); printv(a3);
	printf("b1 = "); printv(b1);
	printf("b2 = "); printv(b2);
	printf("b3 = "); printv(b3);
}

/***************************************/

/* crystal */

void Fcc(){
	/* face-centered cubic */
	int n=0;
	S="fcc";
	N=4*Nx*Ny*Nz;
	R=vec1d(N);
	for(int i=0; i<Nx; i++)
		for(int j=0; j<Ny; j++)
			for(int k=0; k<Nz; k++){
				R[n]=cart(i*a,j*a,k*a);
				R[n+1]=cart(i*a,(j+0.5)*a,(k+0.5)*a);
				R[n+2]=cart((i+0.5)*a,j*a,(k+0.5)*a);
				R[n+3]=cart((i+0.5)*a,(j+0.5)*a,k*a);
				n+=4;
			}
}

void Dia(){
	/* diamond cubic */
	int n=0;
	S="dia";
	N=8*Nx*Ny*Nz;
	R=vec1d(N);
	for(int i=0; i<Nx; i++)
		for(int j=0; j<Ny; j++)
			for(int k=0; k<Nz; k++){
				R[n]=cart(i*a,j*a,k*a);
				R[n+1]=cart(i*a,(j+0.5)*a,(k+0.5)*a);
				R[n+2]=cart((i+0.5)*a,j*a,(k+0.5)*a);
				R[n+3]=cart((i+0.5)*a,(j+0.5)*a,k*a);
				R[n+4]=cart((i+0.25)*a,(j+0.25)*a,(k+0.25)*a);
				R[n+5]=cart((i+0.25)*a,(j+0.75)*a,(k+0.75)*a);
				R[n+6]=cart((i+0.75)*a,(j+0.25)*a,(k+0.75)*a);
				R[n+7]=cart((i+0.75)*a,(j+0.75)*a,(k+0.25)*a);
				n+=8;
			}
}

/***************************************/

/* nbor list */

vec c_to_f(vec c){
	/* Cartesian to fractional */
	vec f;
	f=v(3);
	f.v[0]=mod0(dot(c,b1),1);
	f.v[1]=mod0(dot(c,b2),1);
	f.v[2]=mod0(dot(c,b3),1);
	return f;
}

vec f_to_c(vec f){
	/* fractional to Cartesian */
	vec c=v(3);
	for(int i=0; i<3; i++)
		c.v[i]=f.v[0]*a1.v[i]+f.v[1]*a2.v[i]+f.v[2]*a3.v[i];
	return c;
}

vec Vpbc(int i, int j){
	/* pbc displacement: i->j */
	return f_to_c(c_to_f(minus(R[j],R[i])));
}

void pair(int i){
	/* pair displacements & distances for atom i (single) */
	int j;
	vec v;
	double d;
	for(int n=0; n<nb; n++){
		j=Nb[i][n];
		v=Vpbc(i,j);
		d=mag(v);
		V[i][j]=v; // atom i
		V[j][i]=vmul(v,-1); // nbor j
		D[i][j]=d;
		D[j][i]=d;
	}
}

void repair(){
	/* pair displacements & distances for all atoms (whole) */
	int j;
	vec v;
	double d;
	for(int i=0; i<N; i++)
		for(int n=0; n<nb; n++){
			j=Nb[i][n];
			v=Vpbc(i,j);
			d=mag(v);
			V[i][j]=v;
			D[i][j]=d;
		}
}

void make_nbor(){
	/* make nbor list; cutoff=K */
	vec v;
	double d,eps=1e-3;
	for(int i=0; i<N; i++){
		/* get all distances */
		nb=0;
		for(int j=0; j<N; j++){
			v=Vpbc(i,j); // pbc displacement
			d=mag(v); // pbc distance
			V[i][j]=v;
			D[i][j]=d;
			/* pick nbor */
			if(d>0 && d<K+eps){
				Nb[i][nb]=j;
				nb++;
			}
		}
	}
}

/***************************************/

/* potential: LJ,TE */

double V_all(double V(int)){
	/* total potential in cell (whole) */
	/* V: potential func for single atom */
	double v=0;
	for(int i=0; i<N; i++) v+=V(i); // sum indiv potential
	return 0.5*v; // avoid double-counting
}

double V_lj(int i){
	/* Lennard-Jones potential for atom i (single) */
	/* for fcc Ar */
	double eps=0.01,sig=3.40; // params
	double x,v=0; // potential
	for(int n=0; n<nb; n++){
		x=pow(sig/D[i][Nb[i][n]],6);
		v+=4*eps*x*(x-1);
	}
	return v;
}

double f_cut(double r, double R, double L){
	/* Tersoff cutoff func */
	if(r<R-L) return 1;
	if(r<R+L) return 0.5-0.5*sin(M_PI/2*(r-R)/L);
	else return 0;
}

double V_tersoff(int i){
	/* Tersoff potential for atom i (single) */
	/* for dia Si */
	double v=0,v_ij, // potential
	R=2.75,L=0.1,A=3264.7,B=95.373,lam1=3.2394,lam2=1.3258,lam3=lam2,
	n=22.956,m=3,c=4.8381,d=2.0417,cos_t,cos_t0=0,b_ij,beta=0.33675,zeta,gam=1,
	g,d_ij,d_ik; // params

// flow: cos_t -> g -> zeta -> b_ij -> v_ij -> v (total E)
	for(int j=0; j<N; j++){
		if(i==j) continue; // j!=i
		d_ij=D[i][j]; // pair distance
		if(d_ij>R+L) continue; // ignore far atoms
		zeta=0;
		for(int k=0; k<N; k++){
			// cos_t -> g -> zeta
			if(i==k||j==k) continue; // k!=i,j
			d_ik=D[i][k];
			if(d_ik>R+L) continue;
			cos_t=vcos(V[i][j],V[i][k]); // local angle (env-dependence)
			g=gam*(1+c*c/(d*d)-c*c/(d*d+pow(cos_t-cos_t0,2)));
			zeta+=f_cut(d_ik,R,L)*g*exp(pow(lam3*(d_ij-d_ik),m));
		}
		b_ij=pow(1+pow(beta*zeta,n),-0.5/n);
		v_ij=f_cut(d_ij,R,L)*(A*exp(-lam1*d_ij)-b_ij*B*exp(-lam2*d_ij));
		v+=v_ij;
	}
	return v;
}

/***************************************/

/* optimisation: SD,CG */

void perturb(double eps){
	/* randomly perturb all atoms */
	for(int i=0; i<N; i++) R[i]=add(R[i],randv(3,eps));
	repair();
}

void randomise(double d){
	/* randomise whole lattice within cell (not used) */
	/* d: min sep between atoms */
	int n;
	vec v;
	for(int i=0; i<N; i++){
		do{
			R[i]=cart(uniform(0,Nx*a),uniform(0,Ny*a),uniform(0,Nz*a));
			n=0;
			for(int j=0; j<i; j++)
				if(mag(Vpbc(i,j))<d){n++; break;}
		}while(n>0);
	}
	repair();
}

vec F(int i, double V()){
	/* force on atom i under potential V */
	/* V: potential func for single atom */
	/* math: F=-∇V */
	double eps=1e-6,V0=V(i),dVx,dVy,dVz;
	vec Vi=copy(R[i]);
	R[i]=add(Vi,cart(eps,0,0)); // along x-axis
	pair(i); dVx=V(i)-V0;
	R[i]=add(Vi,cart(0,eps,0)); // along y-axis
	pair(i); dVy=V(i)-V0;
	R[i]=add(Vi,cart(0,0,eps)); // along z-axis
	pair(i); dVz=V(i)-V0;
	R[i]=Vi; pair(i); // orig
	return cart(-dVx/eps,-dVy/eps,-dVz/eps); // force vec
}

vec F_lj(int i){
	/* exact LJ force for atom i */
	int j;
	double eps=0.01,sig=3.40; // params
	double r,x;
	vec Fi=z(3);
	for(int n=0; n<nb; n++){
		j=Nb[i][n];
		r=D[i][j];
		x=sig/r;
		Fi=add(Fi,vmul(V[i][j],2*pow(x,14)-pow(x,8)));
	}
	return vmul(Fi,-24*eps/(sig*sig));
}

void SD(double V()){
	/*** steepest descent under potential V ***/
	/* allow atoms to descent under forces */
	int    n=0; // counter
	vec    g0[N],g1[N],R0[N];
	double eps=1e-6,x,y,b,V0=V_all(V),V1,err=1;
	FILE   *f=fopen("log.xyz","w");

	line();
	while(err>1e-8 && n<100){
		/* calc steepest step size */
		x=0; y=0;
		for(int i=0; i<N; i++){
			R0[i]=copy(R[i]); // orig coord
			g0[i]=F(i,V); // force
			x+=dot(g0[i],g0[i]);
		}
		for(int i=0; i<N; i++) R[i]=add(R0[i],vmul(g0[i],eps));
		repair();

		for(int i=0; i<N; i++){
			g1[i]=F(i,V);
			y+=dot(g1[i],g0[i]);
		}
		b=-eps*x/(y-x); // step size

		/* update new coord */
		for(int i=0; i<N; i++)
			R[i]=add(R0[i],vmul(g0[i],b));
		repair();

		/* error */
		V1=V_all(V); // updated potential
		err=fabs((V1-V0)/V0);
		V0=V1;

		/* coord & potential log */
		printf(" sd %2d: %lf\n",n,V1);
		print_latt(f,"");
		n++;
	}
	line();
	fclose(f);
}

void CG(double V()){
	/*** conjugate gradient under potential V ***/
	int    n=0; // counter
	vec    g0[N],g1[N],g2[N],h[N],R0[N];
	double eps=1e-6,x,y,u,v,b,g,V0=V_all(V),V1,err=1;
	FILE   *f=fopen("log.xyz","w");

	u=0;
	for(int i=0; i<N; i++){
		g2[i]=F(i,V);
		h[i]=copy(g2[i]);
		u+=dot(g2[i],g2[i]);
	}

	line();
	while(err>1e-8 && n<100){
		/* calc step size */
		x=0; y=0;
		for(int i=0; i<N; i++){
			R0[i]=copy(R[i]); // orig coord
			g0[i]=F(i,V); // force
			x+=dot(g0[i],h[i]);
		}
		for(int i=0; i<N; i++) R[i]=add(R0[i],vmul(h[i],eps));
		repair();

		for(int i=0; i<N; i++){
			g1[i]=F(i,V);
			y+=dot(g1[i],h[i]);
		}
		b=-eps*x/(y-x); // step size

		/* update new coord */
		for(int i=0; i<N; i++)
			R[i]=add(R0[i],vmul(h[i],b));
		repair();

		/*********************/

		v=0;
		for(int i=0; i<N; i++){
			g2[i]=F(i,V);
			v+=dot(g2[i],g2[i]);
		}
		g=v/u;
		u=v;

		for(int i=0; i<N; i++) h[i]=add(g2[i],vmul(h[i],g));

		/* error */
		V1=V_all(V); // updated potential
		err=fabs((V1-V0)/V0);
		V0=V1;

		/* coord & potential log */
		printf(" cg %2d: %lf\n",n,V1);
		print_latt(f,"");
		n++;
	}
	line();
	fclose(f);
}

/***************************************/

int main(){
	// [##]
	//Fcc(); // Ar, LJ
	Dia(); // Si, TE
	info();
	init();
	make_nbor();
	/****************************/
/*
	// Ar, LJ [##]
	printf("true: %f\n",V_all(V_lj));
	perturb(1.2);
	printf("init: %f\n",V_all(V_lj));
	//SD(V_lj);
	CG(V_lj);
	printf("conv: %f\n",V_all(V_lj));
*/
	/****************************/
/*
	// Si, TE [##]
	printf("true: %f\n",V_all(V_tersoff));
	perturb(0.2);
	printf("init: %f\n",V_all(V_tersoff));
	//SD(V_tersoff);
	CG(V_tersoff);
	printf("conv: %f\n",V_all(V_tersoff));
*/
	/****************************/
	return 0;
}
