//2D POISSON EQUATION WITH JUMP (INTERIOR) BOUNDARY CONDITIONS. GHOST FLUID METHOD (GFM) IS IMPLEMENTED.
//Ref. - Liu et al. (2000, JCP)
//Example 5 (Fig. 6)
//CELL CENTERED FINITE DIFFERENCE DISCRETIZATION USING GHOST CELLS.
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>
#define TOL 1e-6	//convergence criteria
#define EPS 1e-16
#define SMALL 1e-8
using namespace std;
const double WIDTH0=-1.0,WIDTH1=1.0;	//x range
const double HEIGHT0=-1.0,HEIGHT1=1.0;	//y range
const int I=64,J=64;	//mesh
#include "cfd_solvers.h"
#include "GMG2.cpp"
//#include "MG_FDQGMRES.cpp"
#include "MG_BICGSTAB.cpp"
//class GFM:public MG_FDQGMRES
class GFM:public MG_BICGSTAB
{
	double *Xm,*Ym;	//the mesh variables
	double *CX,*CY;	//the node variables
	double **TH,**f;	//solution field and RHS vector
	double **Phi;	//level set field

	void grid_gen();	//grid generation
	void LS();	//create the level set field
	void f_jump(int f_flag,int i0,int j0,int i1,int j1,double Beta_p,double Beta_m,double *Beta_f,double *F_f);	//calculation of jump conditions at the interior boundary
	void Poisson(int *C,int *R,double *A,double *b);	//create the discretized equations
	public:
			GFM(); ~GFM();	//memory management
			void solve();	//solver
			void write();	//file output for Tecplot360
};
//GFM::GFM():MG_FDQGMRES(8,I*J,1000)	//MG_FDQGMRES parameterized constructor
GFM::GFM():MG_BICGSTAB(I*J,1000)	//MG_BICGSTAB parameterized constructor
{
	Xm=new double[I+1];
	Ym=new double[J+1];
	CX=new double[I+1];
	CY=new double[J+1];
	TH=new double*[J+1];
	f=new double*[J+1];
	Phi=new double*[J+2];
	for(int i=0;i<J+2;i++)
	{
		Phi[i]=new double[I+2];
		if(i<J+1)
		{
			TH[i]=new double[I+1];
			f[i]=new double[I+1];
		}
	}
	cout<<"GFM: MEMORY ALLOCATED"<<endl;
}
GFM::~GFM()
{
	for(int i=0;i<J+2;i++)
	{
		if(i<J+1)
		{
			delete[] TH[i];
			delete[] f[i];
		}
		delete[] Phi[i];
	}
	delete[] Xm;
	delete[] Ym;
	delete[] CX;
	delete[] CY;
	delete[] TH;
	delete[] f;
	delete[] Phi;
	cout<<"GFM: MEMORY RELEASED"<<endl;
}
void GFM::grid_gen()
{
	double dx=(WIDTH1-WIDTH0)/I,dy=(HEIGHT1-HEIGHT0)/J;	//grid size
	Xm[0]=WIDTH0; Ym[0]=HEIGHT0;
	CX[0]=Xm[0]-0.5*dx; CY[0]=Ym[0]-0.5*dy;
	for(int i=1;i<=I;i++)
	{
		Xm[i]=Xm[i-1]+dx;
		CX[i]=CX[i-1]+dx;
	}
	for(int j=1;j<=J;j++)
	{
		Ym[j]=Ym[j-1]+dy;
		CY[j]=CY[j-1]+dy;
	}
}
void GFM::LS()
{
	double Xc=0.0,Yc=0.0,rad=0.5,Dc;	//circle parameters
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			Dc=sqrt(pow((CX[i]-Xc),2.0)+pow((CY[j]-Yc),2.0));	//calculate distance
			Phi[j][i]=Dc-rad;	//calculate level set (outward pointing normal)
		}
	}
}
void GFM::f_jump(int f_flag,int i0,int j0,int i1,int j1,double Beta_p,double Beta_m,double *Beta_f,double *F_f)
{
	double Phi_p,Phi_m;	//positive and negative level sets
	double jump_0,jump_1,jump,norm_0,norm_1;	//jump conditions and component of normal vector
	double dn;	//grid size
	Phi_p=(Phi[j0][i0]>0.0)?Phi[j0][i0]:Phi[j1][i1];
	Phi_m=(Phi[j0][i0]<0.0)?Phi[j0][i0]:Phi[j1][i1];
	*Beta_f=Beta_p*Beta_m*(abs(Phi_p)+abs(Phi_m))/(Beta_p*abs(Phi_m)+Beta_m*abs(Phi_p));
	Phi_m=abs(Phi[j1][i1])/(abs(Phi[j0][i0])+abs(Phi[j1][i1]));	//theta (variable is reused)
	jump_0=0.0;	//i0,j0
	jump_1=0.0;	//i1,j1
	jump=jump_0*Phi_m+jump_1*(1.0-Phi_m);	//jump in magnitude
	jump_0=2.0;
	jump_1=2.0;
	if(i0==i1)	//y direction
	{
		dn=abs(Ym[j1]-Ym[j0]);
		norm_0=2.0*CY[j0];	//i0,j0
		norm_1=2.0*CY[j1];	//i1,j1
	}
	else	//x direction
	{
		dn=abs(Xm[i1]-Xm[i0]);
		norm_0=2.0*CX[i0];	//i0,j0
		norm_1=2.0*CX[i1];	//i1,j1
	}
	Phi_p=jump_0*norm_0*Phi_m+jump_1*norm_1*(1.0-Phi_m);	//jump in slope (variable is reused)
	if((Phi[j0][i0]<=0.0)&&(Phi[j1][i1]>0.0)) *F_f=*Beta_f*(jump/(dn*dn)+f_flag*Phi_p*Phi_m/(Beta_p*dn));
	else *F_f=-*Beta_f*(jump/(dn*dn)+f_flag*Phi_p*Phi_m/(Beta_m*dn));
}
void GFM::Poisson(int *C,int *R,double *A,double *b)
{
	double dx=Xm[1]-Xm[0],dy=Ym[1]-Ym[0];	//grid size
	double alp=pow((dx/dy),2.0);	//grid aspect ratio
	double Beta_p=1.0,Beta_m=1.0;	//coefficients in the two subdomains
	double TH_l[J+1],TH_r[J+1],TH_t[I+1],TH_b[I+1];	//Dirichlet boundary conditions
	double A_e,A_w,A_n,A_s,F_e,F_w,F_n,F_s;	//coefficients and jump conditions of the Poisson equation
	for(int j=1;j<=J;j++)	//left and right boundary conditions
		TH_l[j]=TH_r[j]=1.0+log(2.0*sqrt(1.0+CY[j]*CY[j]));
	for(int i=1;i<=I;i++)	//bottom and top boundary conditions
		TH_b[i]=TH_t[i]=1.0+log(2.0*sqrt(1.0+CX[i]*CX[i]));
	int cnt=0; R[cnt]=0;
	for(int j=1;j<=J;j++)	//create source term
		for(int i=1;i<=I;i++)
			if(Phi[j][i]<0.0) f[j][i]=0.0;
			else f[j][i]=0.0;
	for(int j=1;j<=J;j++)	//create discreized linear system
	{
		for(int i=1;i<=I;i++)
		{
			if(SGN(Phi[j][i]*Phi[j][i-1])==1)	//left cell face
			{
				A_w=(Phi[j][i]>0.0)?Beta_p:Beta_m;
				F_w=0.0;
			}
			else f_jump(-1,i,j,i-1,j,Beta_p,Beta_m,&A_w,&F_w);
			if(SGN(Phi[j][i]*Phi[j][i+1])==1)	//right cell face
			{
				A_e=(Phi[j][i]>0.0)?Beta_p:Beta_m;
				F_e=0.0;
			}
			else f_jump(1,i,j,i+1,j,Beta_p,Beta_m,&A_e,&F_e);
			if(SGN(Phi[j][i]*Phi[j-1][i])==1)	//bottom cell face
			{
				A_s=(Phi[j][i]>0.0)?Beta_p:Beta_m;
				F_s=0.0;
			}
			else f_jump(-1,i,j,i,j-1,Beta_p,Beta_m,&A_s,&F_s);
			if(SGN(Phi[j][i]*Phi[j+1][i])==1)	//top cell face
			{
				A_n=(Phi[j][i]>0.0)?Beta_p:Beta_m;
				F_n=0.0;
			}
			else f_jump(1,i,j,i,j+1,Beta_p,Beta_m,&A_n,&F_n);
			b[(j-1)*I+i-1]=dx*dx*(f[j][i]+F_e+F_w+F_n+F_s);
			A_n*=alp; A_s*=alp;
			if((i==1)&&(j==1))	//bottom left corner cell
			{
				b[(j-1)*I+i-1]-=2.0*(A_w*TH_l[1]+A_s*TH_b[1]);
				C[cnt]=(j-1)*I+i-1;	//1,1 term
				A[cnt]=-(A_e+A_n+2.0*A_w+2.0*A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,1 term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==1)&&(j>1)&&(j<J))	//left boundary
			{
				b[(j-1)*I+i-1]-=2.0*A_w*TH_l[j];
				C[cnt]=(j-2)*I+i-1;	//1,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,j term
				A[cnt]=-(A_e+A_n+A_s+2.0*A_w);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,j term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==1)&&(j==J))	//top left corner cell
			{
				b[(j-1)*I+i-1]-=2.0*(A_w*TH_l[J]+A_n*TH_t[1]);
				C[cnt]=(j-2)*I+i-1;	//1,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,J term
				A[cnt]=-(A_e+A_s+2.0*A_n+2.0*A_w);
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,J term
				A[cnt]=A_e;
				cnt++;
			}
			else if((j==J)&&(i>1)&&(i<I))	//top boundary
			{
				b[(j-1)*I+i-1]-=2.0*A_n*TH_t[i];
				C[cnt]=(j-2)*I+i-1;	//i,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,J term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,J term
				A[cnt]=-(A_e+A_w+A_s+2.0*A_n);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,J term
				A[cnt]=A_e;
				cnt++;
			}
			else if((j==J)&&(i==I))	//top right corner cell
			{
				b[(j-1)*I+i-1]-=2.0*(A_e*TH_r[J]+A_n*TH_t[I]);
				C[cnt]=(j-2)*I+i-1;	//I,J-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,J term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,J term
				A[cnt]=-(A_w+A_s+2.0*A_e+2.0*A_n);
				cnt++;
			}
			else if((i==I)&&(j>1)&&(j<J))	//right boundary
			{
				b[(j-1)*I+i-1]-=2.0*A_e*TH_r[j];
				C[cnt]=(j-2)*I+i-1;	//I,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,j term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,j term
				A[cnt]=-(A_w+A_n+A_s+2.0*A_e);
				cnt++;
				C[cnt]=j*I+i-1;	//I,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((i==I)&&(j==1))	//bottom right corner cell
			{
				b[(j-1)*I+i-1]-=2.0*(A_e*TH_r[1]+A_s*TH_b[I]);
				C[cnt]=(j-1)*I+i-2;	//I-1,1 term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,1 term
				A[cnt]=-(A_w+A_n+2.0*A_e+2.0*A_s);
				cnt++;
				C[cnt]=j*I+i-1;	//I,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else if((j==1)&&(i>1)&&(i<I))	//bottom boundary
			{
				b[(j-1)*I+i-1]-=2.0*A_s*TH_b[i];
				C[cnt]=(j-1)*I+i-2;	//i-1,1 term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,1 term
				A[cnt]=-(A_e+A_w+A_n+2.0*A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,1 term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,2 term
				A[cnt]=A_n;
				cnt++;
			}
			else	//inner domain
			{
				C[cnt]=(j-2)*I+i-1;	//i,j-1 term
				A[cnt]=A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,j term
				A[cnt]=A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,j term
				A[cnt]=-(A_e+A_w+A_n+A_s);
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,j term
				A[cnt]=A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,j+1 term
				A[cnt]=A_n;
				cnt++;
			}
			R[(j-1)*I+i]=cnt;
		}
	}
}
void GFM::solve()
{
	grid_gen();
	LS();
	Poisson(C0,R0,A0,b0);	//create the discrete equations
	GMG2::setup(I/2,J/2,C0,R0,A0,C1,R1,A1);	//create the coarse grid coefficient matrices
	GMG2::setup(I/4,J/4,C1,R1,A1,C2,R2,A2);
	GMG2::setup(I/8,J/8,C2,R2,A2,C3,R3,A3);
	GMG2::setup(I/16,J/16,C3,R3,A3,C4,R4,A4);
	//MG_FDQGMRES::solve(C0,R0,A0,X0,b0);
	MG_BICGSTAB::solve(C0,R0,A0,X0,b0);
	for(int j=1;j<=J;j++)	//update the solution matrix
		for(int i=1;i<=I;i++)
			TH[j][i]=X0[IND(I,i,j)];
}
void GFM::write()
{
	ofstream p_out("TH_5.dat");
	p_out<<"TITLE = \"Solution field\""<<endl;
	p_out<<"VARIABLES = \"CX\",\"CY\",\"TH\""<<endl;
	p_out<<"ZONE I="<<I<<", J="<<J<<", DATAPACKING=BLOCK"<<endl;
	for(int j=1;j<=J;j++)	//print X co-ordinates of cell centers
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<CX[i];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=1;j<=J;j++)	//print Y co-ordinates of cell centers
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<CY[j];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<TH[j][i];
		p_out<<endl;
	}
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
