//1D POISSON EQUATION WITH JUMP (INTERIOR) BOUNDARY CONDITIONS. GHOST FLUID METHOD (GFM) IS IMPLEMENTED.
//Ref. - Liu et al. (2000, JCP)
//Example 1 (1st part, Fig. 1)
//CELL CENTERED FINITE DIFFERENCE DISCRETIZATION USING GHOST CELLS.
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>
#define TOL 1e-6	//convergence criteria
#define EPS 1e-16
using namespace std;
const double WIDTH=1.0;	//domain size
const int I=64;	//mesh
#include "FDQGMRES_1d.cpp"
#include "cfd_solvers.h"
class GFM:public FDQGMRES
{
	double TH[I+1],f[I+1];	//solution field and RHS vector
	double Phi[I+2];	//level set field
	void LS();	//create the level set field
	void Poisson(int *C,int *R,double *A,double *b);	//create the discretized equations
	public:
			GFM();
			void solve();	//solver
			void write();	//file output for Tecplot360
};
GFM::GFM():FDQGMRES(8,I*I,1000)	//FDQGMRES parameterized constructor
{
}
void GFM::LS()
{
	double dx=WIDTH/I;	//grid size
	for(int i=1;i<=I;i++) Phi[i]=(i-0.5)*dx-0.5;
	Phi[0]=Phi[1]; Phi[I+1]=Phi[I];
}
void GFM::Poisson(int *C,int *R,double *A,double *b)
{
	double dx=WIDTH/I;	//grid size
	double Phi_p,Phi_m;	//positive and negative level sets
	double Beta_p=1.0,Beta_m=1.0;	//coefficients
	double a_j=1.0,b_j=0.0;	//jump conditions
	double TH_l=0.0,TH_r=2.0;	//Dirichlet boundary conditions
	double A_e,A_w,F_l,F_r;
	int cnt=0; R[cnt]=0;
	for(int i=1;i<=I;i++)	//create source term
		f[i]=0.0;
	for(int i=1;i<=I;i++)	//create discretized linear system
	{
		if(SGN(Phi[i]*Phi[i+1])==1)
		{
			A_e=(Phi[i]>0.0)?Beta_p:Beta_m;
			F_r=0.0;
		}
		else
		{
			Phi_p=(Phi[i]>0.0)?Phi[i]:Phi[i+1];
			Phi_m=(Phi[i]<0.0)?Phi[i]:Phi[i+1];
			A_e=Beta_p*Beta_m*(abs(Phi_p)+abs(Phi_m))/(Beta_p*abs(Phi_m)+Beta_m*abs(Phi_p));
			Phi_p=abs(Phi[i+1])/(abs(Phi[i])+abs(Phi[i+1]));	//theta (variable is reused)
			if((Phi[i]<=0.0)&&(Phi[i+1]>0.0)) F_r=A_e*(a_j+b_j*Phi_p*dx/Beta_p);
			else F_r=-A_e*(a_j+b_j*Phi_p*dx/Beta_m);
		}
		if(SGN(Phi[i]*Phi[i-1])==1)
		{
			A_w=(Phi[i]>0.0)?Beta_p:Beta_m;
			F_l=0.0;
		}
		else
		{
			Phi_p=(Phi[i]>0.0)?Phi[i]:Phi[i-1];
			Phi_m=(Phi[i]<0.0)?Phi[i]:Phi[i-1];
			A_w=Beta_p*Beta_m*(abs(Phi_p)+abs(Phi_m))/(Beta_p*abs(Phi_m)+Beta_m*abs(Phi_p));
			Phi_m=abs(Phi[i-1])/(abs(Phi[i])+abs(Phi[i-1]));	//theta (variable is reused)
			if((Phi[i]<=0.0)&&(Phi[i-1]>0.0)) F_l=A_w*(a_j-b_j*Phi_m*dx/Beta_p);
			else F_l=-A_w*(a_j-b_j*Phi_m*dx/Beta_m);
		}
		b[i-1]=dx*dx*f[i]+F_l+F_r;
		if(i==1)	//left boundary
		{
			b[i-1]-=2.0*A_w*TH_l;
			C[cnt]=i-1;	//1 term
			A[cnt]=-(A_e+2.0*A_w);
			cnt++;
			C[cnt]=i;	//2 term
			A[cnt]=A_e;
			cnt++;
		}
		else if(i==I)	//right boundary
		{
			b[i-1]-=2.0*A_e*TH_r;
			C[cnt]=i-2;	//I-1 term
			A[cnt]=A_w;
			cnt++;
			C[cnt]=i-1;	//I term
			A[cnt]=-(A_w+2.0*A_e);
			cnt++;
		}
		else	//inner domain
		{
			C[cnt]=i-2;	//i-1 term
			A[cnt]=A_w;
			cnt++;
			C[cnt]=i-1;	//i term
			A[cnt]=-(A_e+A_w);
			cnt++;
			C[cnt]=i;	//i+1 term
			A[cnt]=A_e;
			cnt++;
		}
		R[i]=cnt;
	}
}
void GFM::solve()
{
	LS();
	Poisson(C0,R0,A0,b0);	//create the discrete equations
	FDQGMRES::solve(C0,R0,A0,X0,b0);
	for(int i=1;i<=I;i++)	//update the solution matrix
		TH[i]=X0[i-1];
}
void GFM::write()
{
	double dx=WIDTH/I;	//grid size
	ofstream p_out("TH_1a.dat");
	p_out<<"0\t0"<<endl;	//left boundary
	for(int i=1;i<=I;i++) p_out<<(i-0.5)*dx<<"\t"<<TH[i]<<endl;	//inner domain
	p_out<<"1.0\t2.0";	//right boundary
	p_out.close();
	cout<<"GFM: WRITE SUCCESSFULL"<<endl;
}
int main()
{
	GFM ms;
	clock_t start = clock();
	ms.solve();
	clock_t end = clock();
	cout<<"Ellapsed time = "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;
	ms.write();
	return 0;
}
