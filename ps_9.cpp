//2D POISSON EQUATION WITH JUMP (INTERIOR) BOUNDARY CONDITIONS. GHOST FLUID METHOD (GFM) IS IMPLEMENTED.
//Ref. - Liu et al. (2000, JCP)
//Example 9 (Fig. 10)
//CELL CENTERED FINITE DIFFERENCE DISCRETIZATION USING GHOST CELLS.
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>
#define TOL 1e-6	//convergence criteria
#define EPS 1e-8
#define SMALL 1e-8
using namespace std;
const double WIDTH0=-1.0,WIDTH1=1.0;	//x range
const double HEIGHT0=0.0,HEIGHT1=3.0;	//y range
const int I=64,J=96;	//mesh
#include "cfd_solvers.h"
#include "NRBS.cpp"
#include "GMG2.cpp"
//#include "MG_FDQGMRES.cpp"
#include "MG_BICGSTAB.cpp"
inline double func(double theta,double alp);	//function
inline double d_func(double theta,double alp);	//function derivative
//class GFM:public MG_FDQGMRES
class GFM:public MG_BICGSTAB
{
	double *Xm,*Ym;	//the mesh variables
	double *CX,*CY;	//the node variables
	double **TH,**f;	//solution field and RHS vector
	double **Phi;	//level set field
	int **tag;	//tags of interfacial cells

	void grid_gen();	//grid generation
	void LS();	//create the level set field
	void reinit();	//LS reinitialization algorithm (Fast Sweep Method, Ref. Zhao (2004, Math. Comp.)
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
	tag=new int*[J+2];
	for(int i=0;i<J+2;i++)
	{
		Phi[i]=new double[I+2];
		tag[i]=new int[I+2];
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
		delete[] tag[i];
	}
	delete[] Xm;
	delete[] Ym;
	delete[] CX;
	delete[] CY;
	delete[] TH;
	delete[] f;
	delete[] Phi;
	delete[] tag;
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
inline double func(double theta,double alp)
{
	return ((0.6*cos(theta)-0.3*cos(3.0*theta))*tan(alp)-(0.7*sin(theta)-0.07*sin(3.0*theta)+0.2*sin(7.0*theta)));
}
inline double d_func(double theta,double alp)
{
	return ((-0.6*sin(theta)+0.9*sin(3.0*theta))*tan(alp)-(0.7*cos(theta)-0.21*cos(3.0*theta)+1.4*cos(7.0*theta)));
}
void GFM::LS()
{
	double Xc,Yc,rad,Dc,alpha,theta_l,theta_r,x,y;
	double F_l,F_alp,F_r;
	int n_pts;	//for ray tracing algorithm
	double op,pc;	//distances-> o: center of the curve, c: cell center, and p: intersection points
	Xc=0.0; Yc=1.5;	//center of curve
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			n_pts=0;	//reinitialization
			Phi[j][i]=10.0;	//initial guess
			Dc=sqrt(pow((CX[i]-Xc),2.0)+pow((CY[j]-Yc),2.0));	//calculate distance
			alpha=atan2((CY[j]-Yc),(CX[i]-Xc));

			theta_r=0.0;	//scan entire domain

			while(theta_r<2.0*M_PI)
			{
				theta_l=theta_r;
				F_l=func(theta_l,alpha);
				do
				{
					theta_r+=0.005;
					F_r=func(theta_r,alpha);
				}
				while(SGN(F_l*F_r)==1.0);

				theta_r=NRBS(200,theta_l,theta_r,alpha,func,d_func);
				if(theta_r>2.0*M_PI) break;
				x=0.6*cos(theta_r)-0.3*cos(3.0*theta_r);
				y=Yc+0.7*sin(theta_r)-0.07*sin(3.0*theta_r)+0.2*sin(7.0*theta_r);
				rad=sqrt(pow((x-Xc),2.0)+pow((y-Yc),2.0));
				F_r=Dc-rad;
				Phi[j][i]=(abs(F_r)<abs(Phi[j][i]))?F_r:Phi[j][i];	//store minimum value (value of LS)

				op=sqrt(pow((x-Xc),2.0)+pow((y-Yc),2.0));	//ray tracing algorithm
				pc=sqrt(pow((x-CX[i]),2.0)+pow((y-CY[j]),2.0));
				if(abs(Dc-op-pc)<=1e-10) n_pts++;

				theta_r+=0.005;
			}
			if(n_pts%2==1) Phi[j][i]=abs(Phi[j][i]);	//sign of LS
			else Phi[j][i]=-abs(Phi[j][i]);
		}
	}
}
void GFM::reinit()
{
	double h=MAX2((Xm[1]-Xm[0]),(Ym[1]-Ym[0]));
	double temp,a,b,gamma=0.5*h;	//gamma is the distance parameter
//-------------INITIALIZATION SCHEME (BASED ON INITIAL LS FIELD)----------------------------------------
	for(int j=0;j<=J+1;j++)	//reinitialize the tag values
		for(int i=0;i<=I+1;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//tagging in horizontal direction (inner domain)
	{
		for(int i=1;i<=I;i++)
		{
			if(SGN(Phi[j][i]*Phi[j][i+1])!=1)	//LS changes sign
			{
				if(abs(Phi[j][i])<=abs(Phi[j][i+1])) tag[j][i]=1;	//tag cell having minimum LS value
				else if(i!=I) tag[j][i+1]=1;
			}
		}
	}
	for(int i=1;i<=I;i++)	//tagging in vertical direction (inner domain)
	{
		for(int j=1;j<=J;j++)
		{
			if(SGN(Phi[j][i]*Phi[j+1][i])!=1)	//LS changes sign
			{
				if(abs(Phi[j][i])<=abs(Phi[j+1][i])) tag[j][i]=1;	//tag cell having minimum LS value
				else if(j!=J) tag[j+1][i]=1;
			}
		}
	}
	for(int j=0;j<=J+1;j++)	//reset LS in untagged cells (including ghost cells)
		for(int i=0;i<=I+1;i++)
			if(tag[j][i]!=1) Phi[j][i]=SGN(Phi[j][i])*100.0;
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phi[j][0]=Phi[j][1]; tag[j][0]=tag[j][1];
		Phi[j][I+1]=Phi[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phi[0][i]=Phi[1][i]; tag[0][i]=tag[1][i];
		Phi[J+1][i]=Phi[J][i]; tag[J+1][i]=tag[J][i];
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
void GFM::f_jump(int f_flag,int i0,int j0,int i1,int j1,double Beta_p,double Beta_m,double *Beta_f,double *F_f)
{
	double Phi_p,Phi_m;	//positive and negative level sets
	double jump_0,jump_1,jump;	//jump conditions and component of normal vector
	double dx=Xm[1]-Xm[0],dy=Ym[1]-Ym[0],dn;	//grid size
	VECTOR N0,N1;	//interface normal vectors
	Phi_p=(Phi[j0][i0]>0.0)?Phi[j0][i0]:Phi[j1][i1];
	Phi_m=(Phi[j0][i0]<0.0)?Phi[j0][i0]:Phi[j1][i1];
	*Beta_f=Beta_p*Beta_m*(abs(Phi_p)+abs(Phi_m))/(Beta_p*abs(Phi_m)+Beta_m*abs(Phi_p));
	Phi_m=abs(Phi[j1][i1])/(abs(Phi[j0][i0])+abs(Phi[j1][i1]));	//theta (variable is reused)
	jump_0=-(CX[i0]*CX[i0]+CY[j0]*CY[j0])-exp(CX[i0])*(CX[i0]*CX[i0]*sin(CY[j0])+CY[j0]*CY[j0]);	//i0,j0
	jump_1=-(CX[i1]*CX[i1]+CY[j1]*CY[j1])-exp(CX[i1])*(CX[i1]*CX[i1]*sin(CY[j1])+CY[j1]*CY[j1]);	//i1,j1
	jump=jump_0*Phi_m+jump_1*(1.0-Phi_m);	//jump in magnitude
	N0=VECTOR((0.5*(Phi[j0][i0+1]-Phi[j0][i0-1])/dx),(0.5*(Phi[j0+1][i0]-Phi[j0-1][i0])/dy));
	N0=N0.unit();
	N1=VECTOR((0.5*(Phi[j1][i1+1]-Phi[j1][i1-1])/dx),(0.5*(Phi[j1+1][i1]-Phi[j1-1][i1])/dy));
	N1=N1.unit();
	jump_0=(-20.0*CX[i0]-exp(CX[i0])*((CX[i0]*CX[i0]+2.0*CX[i0])*sin(CY[j0])+CY[j0]*CY[j0]))*N0.x
			+(-20.0*CY[j0]-exp(CX[i0])*(CX[i0]*CX[i0]*cos(CY[j0])+2.0*CY[j0]))*N0.y;	//i0,j0
	jump_1=(-20.0*CX[i1]-exp(CX[i1])*((CX[i1]*CX[i1]+2.0*CX[i1])*sin(CY[j1])+CY[j1]*CY[j1]))*N1.x
			+(-20.0*CY[j1]-exp(CX[i1])*(CX[i1]*CX[i1]*cos(CY[j1])+2.0*CY[j1]))*N1.y;	//i1,j1
	if(i0==i1)	//y direction
	{
		Phi_p=jump_0*N0.y*Phi_m+jump_1*N1.y*(1.0-Phi_m);	//jump in slope (variable is reused)
		dn=dy;
	}
	else	//x direction
	{
		Phi_p=jump_0*N0.x*Phi_m+jump_1*N1.x*(1.0-Phi_m);	//jump in slope (variable is reused)
		dn=dx;
	}
	if((Phi[j0][i0]<=0.0)&&(Phi[j1][i1]>0.0)) *F_f=*Beta_f*(jump/(dn*dn)+f_flag*Phi_p*Phi_m/(Beta_p*dn));
	else *F_f=-*Beta_f*(jump/(dn*dn)+f_flag*Phi_p*Phi_m/(Beta_m*dn));
}
void GFM::Poisson(int *C,int *R,double *A,double *b)
{
	double dx=Xm[1]-Xm[0],dy=Ym[1]-Ym[0];	//grid size
	double alp=pow((dx/dy),2.0);	//grid aspect ratio
	double Beta_p=10.0,Beta_m=1.0;	//coefficients in the two subdomains
	double TH_l[J+1],TH_r[J+1],TH_t[I+1],TH_b[I+1];	//Dirichlet boundary conditions
	double A_e,A_w,A_n,A_s,F_e,F_w,F_n,F_s;	//coefficients and jump conditions of the Poisson equation
	for(int j=1;j<=J;j++)	//left and right boundary conditions
	{
		TH_l[j]=-(CX[1]*CX[1]+CY[j]*CY[j]);
		TH_r[j]=-(CX[I]*CX[I]+CY[j]*CY[j]);
	}
	for(int i=1;i<=I;i++)	//bottom and top boundary conditions
	{
		TH_b[i]=-(CX[i]*CX[i]+CY[1]*CY[1]);
		TH_t[i]=-(CX[i]*CX[i]+CY[J]*CY[J]);
	}
	int cnt=0; R[cnt]=0;
	for(int j=1;j<=J;j++)	//create source term
		for(int i=1;i<=I;i++)
			if(Phi[j][i]<0.0) f[j][i]=exp(CX[i])*(2.0+CY[j]*CY[j]+2.0*sin(CY[j])+4.0*CX[i]*sin(CY[j]));
			else f[j][i]=-40.0;
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
	LS(); reinit();	//create the LS field
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
	ofstream p_out("TH_9.dat");
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
