/*
 * 	locpol.c
 * 		Funciones para el polinomio local en C
 * 	NOTES:
 *		- A?adir todos los kernels....
 *		- Mejorar la convoluci?n utilizando simpson...
 *      - Uses linpack.h routines dpoco_, dpodi_, dposl_
 * 	ERRORS:
 */
#include <R.h>
#include <Rmath.h>
#include <R_ext/Linpack.h>


#define CUAD(x) ((x)*(x))
#define TRI(x) ((x)*(x)*(x)) /* Added by Marek Omelka */
#define MAXDEG 10
#define MAXEVALPTS 1500
#define MAXS 100
#define MAXT 10
#define TOLERANCE 1.0e-200
#define NA_VALUE R_NaReal


/*
 * 	kernel.c
 * 		Funciones para el polinomio local en C
 * 	NOTES:
 *		- A?adir todos los kernels....
 *		- Mejorar la convoluci?n utiulizando simpson...
 *		- El selector del kernel deber?a basarse en RK(), as? la
 *		funci?n en R pasa RK() y aqu? se compara con los diferentes
 *		RK() con una cierta tolerancia.
 * 	ERRORS:
 */

/*
 * Some kernels
 */
double gaussK(double x)
/*	gauss kernel	*/
{
	return( exp(-x*x/2)/sqrt(2*M_PI) );
}

double EpaK(double x)
/*	Epanechnikov kernel	*/
{
	if(fabs(x) <= 1) return( 3*(1-x*x)/4 ); else return( 0 );
}

double Epa2K(double x)
/*	Epanechnikov kernel	*/
{
	if(fabs(x) <= 1) return( 3/(4*sqrt(5))*(1 - x*x/5) ); else  return( 0 );
}

double SqK(double x)
/*	Triangle kernel	*/
{
	if(fabs(x) <= 1) return( 1 ); else return( 0 );
}

double TrianK(double x)
/*	Triangle kernel	*/
{
	if(fabs(x) <= 1) return( 1 - fabs(x) ); else return( 0 );
}

double QuartK(double x)
/*	Quartic kernel	*/
{
	if(fabs(x) <= 1) return( 15./16. * CUAD(1 - x*x) ); else return( 0 );
}

double biweigK(double x)
/*	biweight kernel	*/
{
	if(fabs(x) <= 1) return( CUAD(1 - x*x) ); else return( 0 );
}

double TriweigK(double x)
/*  triweight kernel, added by Marek Omelka. */
{
  if(fabs(x) <= 1) return( 35./32.*TRI(1 - x*x) ); else return( 0 );
}

double tricubK(double x)
/*  tricube kernel */
{
  if(fabs(x) <= 1) return( 70./81.*TRI( 1-TRI(fabs(x)) ) ); else return( 0 );
}

double CosK(double x)
/*  cosin kernel */
{
  if(fabs(x) <= 1) return( M_PI/4*cos(M_PI*x/2) ); else return( 0 );
}

/* Selector del kernel */
typedef double (*funPtr)(double Ktype);

funPtr selKernel(int Ktype)
{
	switch(Ktype)
	{
		case 1:
			return( &EpaK );
		case 2:
			return( &Epa2K );
		case 3:
			return( &TrianK );
		case 4:
			return( &QuartK );
		case 5:
			return( &biweigK );
		case 6:
      return( &TriweigK ); /* added by Marek Omelka. */
		case 7:
			return( &tricubK );
		case 9:
 			return( &CosK );
		case 10:
			return( &SqK );
    case 8:
		default:
			return( &gaussK);
	}
}


double Kconvol(funPtr kernel,double x)
/*
 * 	- Uses simple summation to compute kernel convolution...
 * 	kernel = kernel to covolve.
 * 	x = pt. where convol. jhas to be computed.
 */
{
	double minx = -10.0;
	double maxx = 10.0;
	double delta = 0.1;
	double u,res;

	res = 0.;
	u = minx;
	while(u < maxx)
	{
		res += delta * kernel(u) * kernel(x-u);
		u += delta;
	}
	return( res );
}


/*
 * 	lpest.c
 * 		Funciones para el polinomio local en C
 * 	NOTES:
 *		- Por fin sali? el caso general.
 * 	ERRORS:
 */

/*
 *	Est. densidad.
*/
void parzenRossen(	double *xeval, int *neval,
					double *x, double *weig, int*n, double *bw,
					int *Ktype, double *res)
/*	OK
	xeval = eval. points.
	neval = number of  eval. points.
	x = x data vector.
	weig = vector of weigths for observations.
	n = number of data obs. length of 'x'.
	bw = bandwidth.
	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
	res = result vector.
*/
{
	int i,j;
	funPtr kernel;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		res[i] = 0.;
		for(j=0;j<*n;j++)
		{
			res[i] += weig[j] * kernel( (x[j]-xeval[i]) / (*bw) );
		}
		res[i] /= ((*n) * (*bw));
	}
}


/*
 * smoothers
 */
void simpleSmoother(double *xeval, int *neval,
					double *x, double *y, double *weig, int*n, double *bw,
					int *Ktype, double *res)
/*
 * 	OK
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	res = result vector
 */
{
	int i,j;
	funPtr kernel;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		res[i] = 0.;
		for(j=0;j<*n;j++)
		{
			res[i] += weig[j] * y[j] * kernel( (x[j]-xeval[i]) / (*bw) );
		}
		res[i] /= ((*n) * (*bw));
	}
}


/*
 * loc. pol. estimators.
 */
 void locWeightsEvalxx(double *lpweig, int *neval,
					double *y, int*n, double *res)
/*
	lpweig = neval x n matrix
	neval = number of  eval. points.
	y = y data vector.
	den = vector el denominador de cada estimaci?n.
	res = result matrix, len(xeval) x n dimension.
*/
{
	int i,j;
	/*QUITAR:double wsum;*/

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		res[i] = 0.;
		for(j=0;j<*n;j++)
			if( res[j*(*neval)+i] == NA_VALUE )
				res[i] =  NA_VALUE;
			else
				res[i] += lpweig[j*(*neval)+i] * y[j];
	}
}
void locCteSmoother(double *xeval, int *neval,
					double *x, double *y, double *weig, int*n, double *bw,
					int *Ktype, double *den, double *res)
/*
	xeval = eval. points.
	neval = number of  eval. points.
	x = x data vector.
	y = y data vector.
	weig = vector of weigths for observations.
	n = number of data obs. length of 'x' and 'y'.
	bw = bandwidth.
	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
	den = vector el denominador de cada estimaci?n.
	res = result vector.
*/
{
	int i,j;
	funPtr kernel;
	double s0,t0;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		s0 = t0 = 0.;
		for(j=0;j<*n;j++)
		{
			t0 += weig[j] * y[j] * kernel( (x[j]-xeval[i]) / (*bw) );
			s0 += weig[j] * kernel( (x[j]-xeval[i]) / (*bw) );
		}
		den[i] = s0;
		if( fabs(s0)>TOLERANCE )
			res[i] = t0/s0;
		else
			res[i] = NA_VALUE;
	}
}


void locCteWeights(double *xeval, int *neval,
					double *x, double *weig, int*n, double *bw,
					int *Ktype, double *den, double *res)
/*
	xeval = eval. points.
	neval = number of  eval. points.
	x = x data vector.
	weig = vector of weigths for observations.
	n = number of data obs. length of 'x' and 'y'.
	bw = bandwidth.
	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
	den = vector el denominador de cada estimaci?n.
	res = result matrix, len(xeval) x n dimension.
*/
{
	int i,j;
	funPtr kernel;
	double s0sum;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		s0sum = 0.;
		for(j=0;j<*n;j++)
		{
			res[j*(*neval)+i] = weig[j] * kernel( (x[j]-xeval[i]) / (*bw) );
			s0sum += res[j*(*neval)+i];
		}
		den[i] = s0sum;
		if( fabs(s0sum)>TOLERANCE )
			for(j=0;j<*n;j++)
				res[j*(*neval)+i] /= s0sum;
		else
			for(j=0;j<*n;j++)
				res[j*(*neval)+i] = NA_VALUE;
	}
}


void locLinSmoother(double *xeval, int *neval,
					double *x, double *y, double *weig, int*n, double *bw,
					int *Ktype, double *den, double *beta0, double *beta1)
/*
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	den = vector el denominador de cada estimaci?n.
 * 	res = result vector.
 */
{
	int i,j;
	funPtr kernel;
	double s0,s1,s2,t0,t1,auxK,aux;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		s0 = s1 =s2 = t0 = t1 = 0.;
		for(j=0;j<*n;j++)
		{
			aux = (x[j]-xeval[i])/(*bw);
			auxK = weig[j] * kernel( aux );
			t0 += y[j] * auxK;
			t1 += y[j] * auxK * aux;
			s0 += auxK;
			s1 += auxK * aux;
			s2 += auxK * CUAD(aux);
		}
		den[i] = (s2*s0 - s1*s1);
		if( fabs(den[i])> TOLERANCE )
		{
			beta0[i] = (s2*t0 - s1*t1)/den[i];
			beta1[i] = -(s1*t0 - s0*t1)/den[i];
      beta1[i] /= (*bw);
		}
		else
			beta0[i] = beta1[i] = NA_VALUE;
	}
}

void locLinWeights(double *xeval, int *neval,
					double *x, double *weig, int*n, double *bw,
					int *Ktype, double *den, double *res)
/*
	xeval = eval. points.
	neval = number of  eval. points.
	x = x data vector.
	weig = vector of weigths for observations.
	n = number of data obs. length of 'x' and 'y'.
	bw = bandwidth.
	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
	den = vector el denominador de cada estimaci?n.
	res = result matrix, len(xeval) x n dimension.
*/
{
	int i,j;
	funPtr kernel;
	double wsum,s1,s2,aux,auxK;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		wsum = s1 = s2 = 0.;
		for(j=0;j<*n;j++)
		{
			aux = res[j*(*neval)+i] = (x[j]-xeval[i])/(*bw);
			auxK = weig[j] * kernel( aux );
			s1 += auxK * aux;
			s2 += auxK * CUAD(aux);
		}
		for(j=0;j<*n;j++)
		{
			aux = res[j*(*neval)+i];
			auxK = kernel( aux );
			res[j*(*neval)+i] = weig[j] *( s2 * auxK - s1 * auxK * aux );
			wsum += res[j*(*neval)+i];
		}
		den[i] = wsum;
		if( fabs(wsum)>TOLERANCE )
			for(j=0;j<*n;j++)
				res[j*(*neval)+i] /= wsum;
		else
			for(j=0;j<*n;j++)
				res[j*(*neval)+i] = NA_VALUE;
	}
}

void locCuadSmoother(double *xeval, int *neval,
					double *x, double *y, double *weig, int*n, double *bw,
					int *Ktype, double *den,
					double *beta0, double *beta1, double *beta2)
/*
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	beta = results.
 */
{
	int i,j,k,l,deg;
	funPtr kernel;
	double s[MAXS],t[MAXT],aux,Kaux;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	deg = 2;
	for(i=0;i<*neval;i++)
	{
		/* Inicializaci?n s,t */
		for(k=0;k<=deg;k++)
		{
			t[k] = 0.;
			s[k] = 0.;
		}
		for(l=k;l<=2*deg;l++)
			s[l] = 0.;
		/* C?lculo s,t */
		for(j=0;j<*n;j++)
		{
			aux = (x[j]-xeval[i])/(*bw);
			Kaux = weig[j] * kernel( aux );
			for(k=0;k<=deg;k++)
			{
				t[k] += y[j] * Kaux * R_pow_di(aux,k);
				s[k] += Kaux * R_pow_di(aux,k);
			}
			/* deg = 2*deg;*/
			for(l=k;l<=2*deg;l++)
				s[l] += Kaux * R_pow_di(aux,l);
		}
		den[i] = R_pow_di(s[2],3) + s[0]*R_pow_di(s[3],2) +
				 R_pow_di(s[1],2)*s[4] - s[2]*(2*s[1]*s[3] + s[0]*s[4]);
		if( fabs(den[i]) > TOLERANCE )
		{

			beta0[i] = (R_pow_di(s[3],2)*t[0] - s[2]*s[4]*t[0] + s[1]*s[4]*t[1] +
				R_pow_di(s[2],2)*t[2] - s[3]*(s[2]*t[1] + s[1]*t[2])) / den[i];

			beta1[i] = (s[1]*s[4]*t[0] + R_pow_di(s[2],2)*t[1] - s[2]*(s[3]*t[0] +
				s[1]*t[2]) + s[0]*(-(s[4]*t[1]) + s[3]*t[2])  ) / den[i] ;
			beta1[i] /= (*bw);

			beta2[i] = (R_pow_di(s[2],2)*t[0] - s[1]*s[3]*t[0] + s[0]*s[3]*t[1] +
				R_pow_di(s[1],2)*t[2] - s[2]*(s[1]*t[1] + s[0]*t[2])) / den[i];
			beta2[i] /= (*bw)*(*bw)/2;
		}
		else
			beta0[i] = beta1[i] = beta2[i] = NA_VALUE;
		den[i] = -den[i];
	}
}


void lsSolve(	double*a, int*lda, int*n, double*b, int*bcol,
				double*rcond, double*z, int*info)
/*
 *	a = symmetric matrix...
 *	lda =
 *	n =
 *	b =
 *	bcol =
 *	rcond =
 *	z =
 *	info =
 */
{
	int j;

	dpoco_(a, lda, n, rcond, z, info);
	if( fabs(*rcond)<TOLERANCE)
		warning("Bad conditioned matrix.");
	else if((*info)!=0)
		warning("Bad info result.");
	else
		for(j=0;j<(*bcol);j++)
			dposl_(a,lda,n,&b[j*(*bcol)]);
}


void locPolSmoother(double *xeval, int *neval,
					double *x, double *y, double *weig, int*n,
					double *bw, int *deg, int *Ktype, int *compDet,
					double *den, double *beta)
/*
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	den = det of s matrix.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 *	compDet = 1 computes determinant of s matrix(apox. n h f(x)) en den.
 * 	beta = results.
 *	NOTAS:
 *		El puto Fortran guarda las matrices por columnas!!!
 *		Vaya foll?n......
 */
{
	int i, j, k, l, info, sncol, job;
	funPtr kernel;
	double  s[MAXS], t[MAXT], z[MAXDEG], det[2], aux, Kaux, rcond;

	sncol = (*deg)+1;
	//__TODO__ remove tncol = sncol;
	/*	Selector de kernel */
	kernel = selKernel(*Ktype);
	/*	operaci?n */
	for(i=0;i<(*neval);i++)
	{
		/* Inicializaci?n s,t */
		for(k=0;k<=(*deg);k++)
		{
			t[k] = 0.;
			for(l=k;l<=(*deg);l++)
					s[l*sncol+k] = 0.;
		}
		/* C?lculo s,t */
		for(j=0;j<(*n);j++)
		{
			aux = (x[j]-xeval[i])/(*bw);
			Kaux = weig[j] * kernel( aux );
			for(k=0;k<=(*deg);k++)
			{
				t[k] += y[j] * Kaux * R_pow_di(aux,k);
				for(l=k;l<=(*deg);l++)
					s[l*sncol+k] += Kaux * R_pow_di(aux,k+l);	/* s[k*sncol+l]=*/
			}
		}
		/* Compute sol.  (det ?) */
		rcond=TOLERANCE+1;
		info=0;
		dpoco_(s, &sncol, &sncol, &rcond, z, &info);
		if( fabs(rcond)<TOLERANCE)
		{
			warning("Bad conditioned mAAAAAAtrix at %f.\n",xeval[i]);
			for(k=0;k<=(*deg);k++)
				beta[k*(*neval)+i] = NA_VALUE;
			den[i] = NA_VALUE;
		}
		else if( info!=0 )
		{
			warning("Bad info result  at %f.\n",xeval[i]);
			for(k=0;k<=(*deg);k++)
				beta[k*(*neval)+i] = NA_VALUE;
			den[i] = NA_VALUE;
		}
		else
		{
			dposl_(s,  &sncol,  &sncol, t);
			job = 10;
			if( (*compDet) )
			{
				dpodi_(s,  &sncol,  &sncol, det, &job);
				if( 1.0 <= det[0] && det[0] <= 10.0 )
					den[i] = det[0] * R_pow (10.0, det[1]);
				else
					den[i] = det[0];
			}
			/* Colocar la soluci?n en beta y reajustar */
			aux = (*bw);
			l = 1;
			beta[i] = t[0];
			for(k=1;k<=(*deg);k++)
			{
				l *= k;
				beta[k*(*neval)+i] = t[k]*l/aux;
				aux *= (*bw);
			}
		}
	}
}




/*
 * Leave One out computations.
 * Estimation of m(x_i) i=1,...n without ith observation (x_i,y_i).
 */
void looLocPolSmoother(	double *x, double *y, double *weig, int*n,
						double *bw, int *deg, int *Ktype, int *compDet,
						double *den, double *beta)
/*
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	beta = results.
 *	NOTAS:
 *		El puto Fortran guarda las matrices por columnas!!!
 *		Vaya foll?n......
 *	RETURN
 *		Leave one out local polinomial estimator. It returns the estimation
 *	of m(x_i) but without considering the i-th observation, in order to avoid
 *	possible biases.
 */
{
	int i, j, k, l, info, sncol, job;  //__TODO__ remove tncol,
	funPtr kernel;
	double  s[MAXS], t[MAXT], z[MAXDEG], det[2], aux, Kaux, rcond;

	sncol = (*deg)+1;
	//__TODO__ remove tncol = sncol;
	/*	Selector de kernel */
	kernel = selKernel(*Ktype);
	/*	operaci?n */
	for(i=0;i<(*n);i++)
	{
		/* Inicializaci?n s,t */
		for(k=0;k<=(*deg);k++)
		{
			t[k] = 0.;
			for(l=k;l<=(*deg);l++)
					s[l*sncol+k] = 0.;
		}
		/* C?lculo s,t */
		for(j=0;j<(*n);j++)
		{
			if( i!=j )	/* Do not consider data x_j to estim. m(x_i) if i==j */
			{
				aux = (x[j]-x[i])/(*bw);
				Kaux = weig[j] * kernel( aux );
				for(k=0;k<=(*deg);k++)
				{
					t[k] += y[j] * Kaux * R_pow_di(aux,k);
					for(l=k;l<=(*deg);l++)
						s[l*sncol+k] += Kaux * R_pow_di(aux,k+l);	/* s[k*sncol+l]=*/
				}
			}
		}
		/* Compute sol.  (det ?) */
		rcond=TOLERANCE+1;
	  info=0;
	  dpoco_(s, &sncol, &sncol, &rcond, z, &info);
		if( fabs(rcond)<TOLERANCE)
		{
			warning("Bad conditioned matrix at %f.\n",x[i]);
			for(k=0;k<=(*deg);k++)
				beta[k*(*n)+i] = NA_VALUE;
			den[i] = NA_VALUE;
		}
		else if( info!=0 )
		{
			warning("Bad info result  at %f.\n",x[i]);
			for(k=0;k<=(*deg);k++)
				beta[k*(*n)+i] = NA_VALUE;
			den[i] = NA_VALUE;
		}
		else
		{
			dposl_(s,  &sncol,  &sncol, t);
			job = 10;
			if( (*compDet) )
			{
				dpodi_(s,  &sncol,  &sncol, det, &job);
				if( 1.0 <= det[0] && det[0] <= 10.0 )
					den[i] = det[0] * R_pow (10.0, det[1]);
				else
					den[i] = det[0];
			}
			/* Colocar la soluci?n en beta y reajustar */
			aux = (*bw);
			l = 1;
			beta[i] = t[0];
			for(k=1;k<=(*deg);k++)
			{
				l *= k;
				beta[k*(*n)+i] = t[k]*l/aux;
				aux *= (*bw);
			}
		}
	}
}



/*
 * 	cv.c
 * 		Funciones para Validaci?n Cruzada para el polinomio local en C.
 * 	NOTES:
 *		- Se proponen evaluadores, que se deben de llamar desde R y
 *		pasarselos a la funci?n 'optim()'.
 * 	ERRORS:
 */

void denCVBwEval(	double *h,double *x, double *weig, int *n,
					int *Ktype,double *res)
/*
 * 	Compute CV evaluation for a given h.
 * 	x = data.
 *	weig = weights.
 *	n = length of data points.
 *	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 *	res = result value.
 *	NOTES:
 *		Cuidado con los pesos, no han de estar normalizados a n.
 */
{
	funPtr kernel;
	int i,j;
	double aux,sumCuadWeig;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	Operation */
	*res = 0.;
	sumCuadWeig = 0.;
	for(i=0;i<((*n)-1);i++)
	{
		for(j=i+1;j<(*n);j++)
		{
			aux = (x[i]-x[j])/(*h);
			*res += weig[i] * weig[j] * ( Kconvol(kernel,aux) -
					2 * (*n) * kernel(aux)/((*n)-1) );
		}
		sumCuadWeig += CUAD(weig[i]);
	}
	*res += sumCuadWeig * Kconvol(kernel,0.) / 2;
	*res = 2 * (*res) / (CUAD((*n)) * (*h));
/*
 * 	*res = 0.;
 * 	for(i=0;i<((*n)-1);i++)
 * 		for(j=i+1;j<(*n);j++)
 * 		{
 * 			aux = (x[i]-x[j])/(*h);
 * 			*res += ( Kconvol(kernel,aux) -
 * 					2 * (*n) * kernel(aux)/((*n)-1) );
 * 		}
 * 	*res += (*n) * Kconvol(kernel,0.) / 2;
 * 	*res = 2 * (*res) / (CUAD((*n)) * (*h))
 */
}

void swap(double *xx,double *yy)
{
	double aux;
	aux = *xx;
	*xx = *yy;
	*yy = aux;
}

void regCVBwEval(	double *h,double *x, double *y, double *weig, int *n,
					int *deg, int *Ktype, double *res)
/*
 * 	Compute CV evaluation for beta0 and a given h,.
 *	h = bandwidth.
 * 	x,y = data.
 *	weig = weights.
 *	n = length of data points.
 *	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 *	res = result value.
 *	NOTE:
 *		The behavior of 'weig' weigths has been not tested but considered.
 */
{
	int i, ult, neval,compDet;
	double xeval,beta[MAXDEG],den[150];

	compDet = 0;	/* Should be zero, det is not needed  */
	neval = 1;
	ult = (*n)-1;
	*res = 0.;
	/* Operation */
	for(i=0;i<(*n);i++)
	{
		/* paso dato i al ultimo sitio */
		swap(&x[i],&x[ult]);
		swap(&y[i],&y[ult]);
		/* computo est. loc. */
		xeval = x[ult];
		locPolSmoother(	&xeval, &neval, x, y, weig, &ult,
						h, deg, Ktype, &compDet, den, beta);
		/* Comp. ress */
		*res += CUAD(y[ult] - beta[0]);
		/* poner cada dato en su sitio */
		swap(&x[i],&x[ult]);
		swap(&y[i],&y[ult]);
	}
	*res /= (*n);
}

void regCVBwEvalB(	double *h,double *x, double *y, double *weig, int *n,
					int *deg, int *Ktype, double *res)
/*
 * 	Compute CV evaluation for beta0 and a given h,.
 *	h = bandwidth.
 * 	x,y = data.
 *	weig = weights.
 *	n = length of data points.
 *	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 *	res = result value.
 *	NOTE:
 *		The behavior of 'weig' weigths has been not tested but considered.
 */
{
	int i,compDet;
	double *beta, *den;

	compDet = 0;
	den = (double *) R_alloc( (*n), sizeof(double) );
	beta = (double *) R_alloc( ((*deg)+1)*(*n), sizeof(double) );
	looLocPolSmoother( x, y, weig, n, h, deg, Ktype, &compDet, den, beta);
	*res = 0.;
	/* Operation */
	for(i=0;i<(*n);i++)
		*res += (CUAD(y[i] - beta[i]) * weig[i]);
	*res /= (*n);
}



/*
 * Loc. pol. Cuadratic estimators
 * (only local constant and local linear)
 * Useful for variance estimation.
 */
void simpleSqSmoother(	double *xeval, int *neval,
						double *x, double *y, int*n, double *bw,
						int *Ktype, double *res)
/*	OK
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	res = result vector.
 * 	RETURN
 * 		Simple Squared Smoother, computes a simple smoother but with a
 * 	data and kernel squared.
 */
{
	int i,j;
	funPtr kernel;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		res[i] = 0.;
		for(j=0;j<*n;j++)
		{
			res[i] +=  CUAD( y[j] * kernel((x[j]-xeval[i])/(*bw)) );
		}
		res[i] /= ((*n) * (*bw));
	}
}

void locCteSqSmoother(	double *xeval, int *neval,
						double *x, double *y, double *weig, int*n, double *bw,
						int *Ktype, double *den, double *res)
/*
 * xeval = eval. points.
 * neval = number of  eval. points.
 * x = x data vector.
 * y = y data vector.
 * weig = vector of weigths for observations.
 * n = number of data obs. length of 'x' and 'y'.
 * bw = bandwidth.
 * Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * den = vector el denominador de cada estimaci?n.
 * res = result vector.
 *	NOTAS:
 *		Cuidado, lo que se ponga en y[i] se pondera por weig[i], con lo que
 *	en caso de datos sesgados no hay que poner los residuos ponderados, sino
 *	s?lo los residuos...
 *	RETURN
 *  	Simple Loc Constant Smoother, computes a simple smoother but with a
 *	data and kernel squared.
 *		It is an estimation of var of estimations of m(x),
 *	i.e.:'var(x)RK(kernel)\over f(x)'. */
{
	int i,j;
	funPtr kernel;
	double s0;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/*	operaci?n */
	for(i=0;i<*neval;i++)
	{
		/* calcular s_0 y res[] */
		s0 = 0.;
		for(j=0;j<*n;j++)
		{
			s0 += weig[j] * CUAD( kernel( (x[j]-xeval[i]) / (*bw) ) );
			res[i] += weig[j] * CUAD( y[j] * kernel( (x[j]-xeval[i]) / (*bw) ) );
		}
		den[i] = s0;
		if( fabs(s0)>TOLERANCE )
			res[i] = res[i]/s0;
		else
			res[i] = NA_VALUE;
	}
}

void locLinSqSmoother(	double *xeval, int *neval,
						double *x, double *y, double *weig, int*n, double *bw,
						int *Ktype, double *den, double *res )
/*
 * 	xeval = eval. points.
 * 	neval = number of  eval. points.
 * 	x = x data vector.
 * 	y = y data vector.
 * 	weig = vector of weigths for observations.
 * 	n = number of data obs. length of 'x' and 'y'.
 * 	bw = bandwidth.
 * 	Ktype = tipo de kernel:0=gaussK,1=EpaK,2=Epa2K,3=TrianK}
 * 	den = vector el denominador de cada estimaci?n, sin elevar al cuadrado.
 * 	res = result vector.
 *	NOTAS:
 *		Cuidado, lo que se ponga en y[i] se pondera por weig[i], con lo que
 *	en caso de datos sesgados no hay que poner los residuos ponderados, sino
 *	s?lo los residuos...
 *	RETURN
 *  	Simple Loc linear Smoother, computes a simple smoother but with a
 * data and kernel squared.
 *		It is an estimation of var of estimations of m(x),
 *	i.e.:'var(x)RK(kernel)\over f(x)'.
 */
{
	int i,j;
	funPtr kernel;
	double s0,s1,s2,auxK,aux;

	/*	Selector de kernel */
	kernel = selKernel(*Ktype);

	/* Operaci?n */
	for(i=0;i<*neval;i++)
	{
		/* calcular los s_i */
		s0 = s1 = s2 = 0.;
		for(j=0;j<*n;j++)
		{
			aux = (x[j]-xeval[i])/(*bw);
			auxK = weig[j] * kernel( aux );
			s0 += auxK;
			s1 += auxK * aux;
			s2 += auxK * CUAD(aux);
		}
		/* calcular los res_i */
		res[i] = 0;
		for(j=0;j<*n;j++)
		{
			aux = (x[j]-xeval[i])/(*bw);
			auxK = weig[j] * kernel( aux );
			res[i] += CUAD( auxK * (s2-aux*s1) * y[j] );
		}
		den[i] = (s2*s0 - s1*s1);
		if( fabs(den[i])> TOLERANCE )
			res[i] /= R_pow_di( den[i],2);/*R_pow( den[i],3./2.);*/
		else
			res[i] = NA_VALUE;
	}
}
