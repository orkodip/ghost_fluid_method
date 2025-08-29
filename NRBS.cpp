//NEWTON-RAPHSON + BISECTION METHOD
//Refer to Press et al., Numerical Recipes in C
double NRBS(int NITER,double xl,double xh,double alpha,double (*func)(double,double),double (*d_func)(double,double))
{
	double f,df,dx,dx_old,fh,fl;
	double temp,rts;
	fl=(*func)(xl,alpha);
	fh=(*func)(xh,alpha);
	if(((fl>0.0)&&(fh>0.0))||((fl<0.0)&&(fh<0.0)))
	{
		cout<<"NRBS(): ERROR IN INITIAL ROOT BRACKETING"<<endl;
		return 0;
	}
	if(abs(fl)<=EPS) return xl;
	else if(abs(fh)<=EPS) return xh;
	if(fl>0.0) { temp=xl; xl=xh; xh=temp; }	//orient the search so that f(xl)<0
	dx=dx_old=abs(xh-xl);	//intial bracket sizes
	rts=0.5*(xl+xh);	//initial guess
	f=(*func)(rts,alpha); df=(*d_func)(rts,alpha);
	for(int i=1;i<=NITER;i++)
	{
		if((((rts-xh)*df-f)*((rts-xl)*df-f)>=0.0)||(abs(2.0*f)>abs(dx_old*df)))
		{	//use BS if NR is out of bounds or interval is not decreasing fast enough
			dx_old=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if(abs(rts-xl)<=EPS) return rts;	//negligible change in root
		}
		else	//use NR
		{
			dx_old=dx;
			dx=f/df;
			temp=rts;
			rts-=dx;
			if(abs(temp-rts)<=EPS) return rts;	//negligible change in root
		}
		if(abs(dx)<=EPS) return rts;	//convergence criteria
		f=(*func)(rts,alpha); df=(*d_func)(rts,alpha);	//evaluate function and its derivative
		if(abs(f)<=EPS) return rts;	//convergence criteria
		else if(f<0.0) xl=rts;	//reduce the bracket size
		else if(f>0.0) xh=rts;
	}
	cout<<"NRBS(): NOT CONVERGED WITHIN ITERATION LIMIT"<<endl;
	return 0;
}
