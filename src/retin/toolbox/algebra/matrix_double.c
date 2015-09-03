/*
Copyright © CNRS 2014. 
Authors: David Picard, Philippe-Henri Gosselin, Romain Negrel, Hedi 
Tabia, Jerome Fellus
Contact: picard@ensea.fr

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/
/**
 * \file matrix.c
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#include "matrix_double.h"
#include "vector_double.h"

#define sgn(x) (((x) >= -1E-14) ? (1) : (-1))
#define SQRT2 (1.414213562373095048801688724209)

void givens_double(double x, double y, double * c, double * s)
{
	double norm = sqrt(x * x + y * y);
  	if (norm == 0.0)
  	{
    		*c = 1.0;
    		*s = 0.0;
  	}
  	else
  	{
    		*c = x / norm;
    		*s = y / norm;
  	}
}

/* Calcule C += A * B; */
void	matrix_CpAB_double (double* pC,double* pA,double* pB,size_t n,size_t p,size_t m)
{
	double r;
	double* p1;
	size_t i,k,j;
	for (j=0;j<m;j++)
	{
		p1 = pA;
		for (k=0;k<p;k++)
		{
			r = *pB++;
			for (i=0;i<n;i++)
			{
				pC[i] += r*p1[i];
			}
			p1 += n;
		}
		pC += n;
	}
}

/* Calcule C += A' * B; */
void	matrix_CpAtB_double (double* pC,double* pA,double* pB,size_t n,size_t p,size_t m)
{
	/*double r;*/
	double* p1;
	size_t i/*,k*/,j;
	for (j=0;j<m;j++)
	{
		p1 = pA;
		for (i=0;i<n;i++)
		{
			/*r = 0;
			for (k=0;k<p;k++)
				r += p1[k]*pB[k];
			*pC++ += r;*/
                        *pC++ += vector_ps_double (p1, pB, p);
			p1 += p;
		}
		pB += p;
	}
}

/* Calcule C += A * B' */
void	matrix_CpABt_double (double* C,double* A,double* B,size_t n,size_t p,size_t m)
{
	double r;
	size_t i,j,k;
	for (k=0;k<p;k++)
	{
		for (j=0;j<m;j++)
		{
			r = B[j+k*m];
			for (i=0;i<n;i++)
				C[i+j*n] += r*A[i+k*n];
		}
	}
}

/* Calcule C += a * a' */
void	matrix_Cpaat_double (double* C,const double* a,size_t n)
{
	size_t j;
	for (j=0;j<n;j++)
		vector_addm_double(C+j*n, a[j], a, n);
}

/* Calcule C += A * B * A' */
void	matrix_CpABAt_double (double* C,double* A,double* B,size_t n,size_t m,double* r)
{
	double a;
	size_t i,k,t,j;
	for (j=0;j<n;j++)
	{
		memset (r,0,m*sizeof(double));
		for (k=0;k<m;k++)
		{
			a = A[j+k*n];
			for (t=0;t<m;t++)
				r[t] += a*B[t+k*m];
		}
		for (t=0;t<m;t++)
		{
			a = r[t];
			for (i=0;i<n;i++)
				C[i+j*n] += A[i+t*n]*a;
		}
	}
}

/* Calcule C += A' * B * A */
void	matrix_CpAtBA_double (double* C,double* A,double* B,size_t n,size_t m,double* r)
{
	double a;
	size_t i,k,t,j;
	for (j=0;j<n;j++)
	{
		memset (r,0,m*sizeof(double));
		for (t=0;t<m;t++)
		{
			a = A[t+j*n];
			for (k=0;k<m;k++)
				r[k] += a*B[k+t*m];
		}
		for (i=0;i<n;i++)
		{
			a = 0;
			for (k=0;k<m;k++)
				a += A[k+i*n]*r[k];
			C[i+j*n] += a;
		}
	}
}

/* Calcule C += A * d * A' */
void	matrix_CpAdAt_double (double* C,double* A,double* d,size_t n,size_t m)
{
	double a;
	size_t i,k,j;
	for (k=0;k<m;k++)
	{
		for (j=0;j<n;j++)
		{
			a = d[k]*A[j+k*n];
			for (i=0;i<n;i++)
			{
				C[i+j*n] += a*A[i+k*n];
			}
		}
	}
}

void	matrix_t_double (double* Y,double* X,size_t n,size_t m)
{
	size_t i,j;
	for (j=0;j<m;j++)
	for (i=0;i<n;i++)
		Y[j+i*m] = X[i+j*n];
}

void	matrix_sum_double (double* s,double* X,size_t n,size_t m)
{
	double	r;
	size_t i,j;
	for (j=0;j<m;j++)
	{
		r = 0;
		for (i=0;i<n;i++)
			r += X[i+j*n];
		s[j] = r;
	}
}

void	matrix_sumt_double (double* s,double* X,size_t n,size_t m)
{
	size_t i,j;
	for (i=0;i<n;i++)
		s[i] = 0;
	for (j=0;j<m;j++)
	{
		for (i=0;i<n;i++)
			s[i] += X[i+j*n];
	}
}


void	matrix_meant_double (double* s,const double* X,size_t n,size_t m)
{
	size_t i,j;
        for (i=0;i<n;i++)
            s[i] = 0;
	for (j=0;j<m;j++)
	{
            for (i=0;i<n;i++)
                s[i] += X[i+j*n];
	}
        for (i=0;i<n;i++)
            s[i] = s[i] / ((double)m);
}



double	matrix_l2_double (const double* pa,const double* pb,size_t n,size_t m)
{
    return sqrt(matrix_l2p2_double(pa,pb,n,m));
}

double	matrix_l2p2_double (const double* v1,const double* v2,size_t n,size_t m)
{
	size_t i;
	double s,x;
	s = 0;
	for (i=0;i<n*m;i++)
	{
		x = v1[i] - v2[i];
		s += x*x;
	}
	return s;
}


#define permuter(a,b) { r = (a); (a) = (b); (b) = r; }

size_t	matrix_ll_double (double* L,double* A,size_t n,size_t m)
{
	size_t i,j,k;
	double r;

	memcpy (L,A,n*m*sizeof(double));

	for (k=0;k<m;k++)
	{
		r = L[k+k*n];
		if (r < 1E-5)
			break;
		r = sqrt(r);
		L[k+k*n] = r;
		for (j=k+1;j<n;j++)
			L[j+k*n] /= r;
		for (j=k+1;j<m;j++)
			for (i=j;i<n;i++)
				L[i+j*n] -= L[i+k*n]*L[j+k*n];

/*		printf ("----\n");
		for (j=0;j<n;j++)
		{
			for (i=0;i<n;i++)
				printf ("%f;",A[j+i*n]-L[j+i*n]);
			printf ("\n");
		}*/
	}

	for (i=1;i<m;i++)
	for (j=0;j<i;j++)
		L[j+i*n] = 0;

	return k;
}

size_t	matrix_qr_double (double* pX,double* pR,size_t n,size_t m)
{
	size_t i,j,k,rank;
	double r,x;
	memset (pR,0,m*m*sizeof(double));
	for (j=0;j<m;j++)
	{
		r = 0;
		for (i=0;i<n;i++)
		{
			x = pX[i+j*n];
			r += x*x;
		}
		if (r < 1E-10)
			break;
		r = sqrt(r);
		pR[j+j*m] = r;
		r = 1/r;
		for (i=0;i<n;i++)
			pX[i+j*n] *= r;
		for (k=j+1;k<m;k++)
		{
			r = 0;
			for (i=0;i<n;i++)
				r += pX[i+j*n]*pX[i+k*n];
			for (i=0;i<n;i++)
				pX[i+k*n] -= pX[i+j*n]*r;
			pR[j+k*m] = r;
		}
	}
        rank = j;
	for (;j<m;j++)
	{
		for (i=0;i<n;i++)
			pX[i+j*n] = 0;
	}
        return rank;
}

void	matrix_qtq_double (double* A,double* d,double* e,size_t n)
{
	size_t i,j,k;
	double rho,beta;
	double* r = (double*)malloc(n*sizeof(double));
	double* s = (double*)malloc(n*sizeof(double));
	memset (r,0,n*sizeof(double));
	memset (s,0,n*sizeof(double));
	for (k=0;k<n-2;k++)
	{
		d[k] = A[k+k*n];
		beta = 0.0;
		for (i=k+1;i<n;i++)
		{
			r[i] = A[k+i*n];
			beta += r[i]*r[i];
		}
		beta = sqrt(beta);
		if (beta < 1E-14)
		{
			r[k+1] = sqrt(2);
		}
		else
		{
			if (r[k+1] > -1E-14)
				rho = 1.0;
			else
				rho = -1.0;
			for (i=k+1;i<n;i++)
				r[i] *= rho/beta;
			r[k+1] += 1.0;
			beta *= -rho;
			rho = sqrt(r[k+1]);
			for (i=k+1;i<n;i++)
				r[i] /= rho;
		}
		e[k] = beta;

		beta = 0.0;
		for (i=k+1;i<n;i++)
		{
			s[i] = 0.0;
			for (j=k+1;j<i;j++)
				s[i] += A[j+i*n]*r[j];
			for (j=i;j<n;j++)
				s[i] += A[i+j*n]*r[j];
			beta += s[i]*r[i];
		}
		for (i=k+1;i<n;i++)
			s[i] -= r[i]*beta/2.0;
		for (j=k+1;j<n;j++)
			for (i=0;i<j+1;i++)
				A[i+j*n] = A[i+j*n] - r[i]*s[j] - s[i]*r[j];
		for (i=k+1;i<n;i++)
			A[i+k*n] = r[i];
	}

	d[n-2] = A[n-2+(n-2)*n];
	e[n-2] = A[n-2+(n-1)*n];
	d[n-1] = A[n-1+(n-1)*n];

	r[n-2] = A[n-2+(n-3)*n];
	r[n-1] = A[n-1+(n-3)*n];
	for (i=0;i<n;i++)
	{
		A[i+(n-2)*n] = 0;
		A[i+(n-1)*n] = 0;
	}
	A[n-2+(n-2)*n] = 1;
	A[n-1+(n-1)*n] = 1;

	for (k=n-3;;k--)
	{
		for (i=k+1;i<n;i++)
			r[i] = A[i+k*n];
		for (i=k+1;i<n;i++)
		{
			rho = 0.0;
			for (j=k+1;j<n;j++)
				rho += r[j]*A[j+i*n];
			for (j=k+1;j<n;j++)
				A[j+i*n] -= rho*r[j];
		}
		for (i=0;i<n;i++)
			A[i+k*n] = 0.0;
		A[k+k*n] = 1;
		if (k == 0)
			break;
	}
	free (r);
	free (s);
}

void	matrix_qtq_t_double (double* A,double* d,double* e,size_t n)
{
	size_t i,j,k;
	double rho,beta;
	double* r = (double*)malloc(n*sizeof(double));
	double* s = (double*)malloc(n*sizeof(double));
	memset (r,0,n*sizeof(double));
	memset (s,0,n*sizeof(double));
	for (k=0;k<n-2;k++)
	{
		d[k] = A[k+k*n];
		beta= 0;
		for (i=k+1;i<n;i++)
		{
			r[i] = A[k+i*n];
			beta = beta + r[i]*r[i];
		}
		beta = sqrt(beta);
		if (beta < 1E-14)
		{
			r[k+1] = sqrt(2);
		}
		else
		{
			if (r[k+1] > -1E-14)
				rho = 1;
			else
				rho = -1;
			for (i=k+1;i<n;i++)
				r[i] *= rho/beta;
			r[k+1] += 1;
			beta *= -rho;
			rho = sqrt(r[k+1]);
			for (i=k+1;i<n;i++)
				r[i] /= rho;
		}
		e[k] = beta;

		beta = 0;
		for (i=k+1;i<n;i++)
		{
			s[i] = 0;
			for (j=k+1;j<i;j++)
				s[i] += A[j+i*n]*r[j];
			for (j=i;j<n;j++)
				s[i] += A[i+j*n]*r[j];
			beta += s[i]*r[i];
		}
		for (i=k+1;i<n;i++)
			s[i] -= r[i]*beta/2;
		for (j=k+1;j<n;j++)
			for (i=0;i<j+1;i++)
				A[i+j*n] = A[i+j*n] - r[i]*s[j] - s[i]*r[j];
	}

	d[n-2] = A[n-2+(n-2)*n];
	e[n-2] = A[n-2+(n-1)*n];
	d[n-1] = A[n-1+(n-1)*n];

	free(s);
	free(r);
}

int	matrix_lu_double (double* M,size_t* idx,double* d,size_t n)
{
	size_t i,imax,j,k;
	double big,dum,sum,temp;
	double* v = (double*)malloc(n*sizeof(double));
	*d = 1.0;
	for (i=0;i<n;i++)
	{
		big = 0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(M[i+j*n])) > big)
				big = temp;
		if (fabs(big) < 1E-10) {
			free (v);
			return 0;
		}
		v[i] = 1.0/big;
	}
	for (j=0;j<n;j++)
	{
		for (i=0;i<j;i++)
		{
			sum = M[i+j*n];
			for (k=0;k<i;k++)
				sum -= M[i+k*n]*M[k+j*n];
			M[i+j*n] = sum;
		}
		big = 0.0;
		imax = 0;
		for (i=j;i<n;i++)
		{
			sum = M[i+j*n];
			for (k=0;k<j;k++)
				sum -= M[i+k*n]*M[k+j*n];
			M[i+j*n] = sum;
			if ((dum=v[i]*fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (k=0;k<n;k++)
			{
				dum = M[imax+k*n];
				M[imax+k*n] = M[j+k*n];
				M[j+k*n] = dum;
			}
			*d = -(*d);
			v[imax] = v[j];
		}
		idx[j] = imax;
		if (M[j+j*n] == 0.0) M[j+j*n] = 1E-10;
		if (j != n)
		{
			dum = 1.0/M[j+j*n];
			for (i=j+1;i<n;i++)
				M[i+j*n] *= dum;
		}
	}
	free (v);
	return 1;
}

void matrix_lub_double (double* b,double* M,size_t *idx,size_t n)
{
	int i,ii=-1,ip,j;
	double sum;
	for (i=0;i<n;i++)
	{
		ip = idx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0)
		{
			for (j=ii;j<i;j++)
				sum -= M[i+j*n]*b[j];
		}
		else if (fabs(sum) > 1E-10)
			ii = i;
		b[i] = sum;
	}
	for (i=n-1;i>=0;i--)
	{
		sum = b[i];
		for (j=i+1;j<n;j++)
			sum -= M[i+j*n]*b[j];
		b[i] = sum/M[i+i*n];
	}
}

int matrix_inv_double (double* Y,double* M,size_t n)
{
	double d, *c;
	size_t i,j;
	size_t* idx = (size_t*)malloc(n*sizeof(size_t));
	if (!matrix_lu_double (M,idx,&d,n)) {
		free (idx);
		return 0;
	}
	c = (double*)malloc(n*sizeof(double));
	for (j=0;j<n;j++)
	{
		for (i=0;i<n;i++)
			c[i] = 0.0;
		c[j] = 1.0;
		matrix_lub_double (c,M,idx,n);
		for (i=0;i<n;i++)
			Y[i+j*n] = c[i];
	}
	free (idx);
	free (c);
	return 1;
}

int matrix_invm_double (double* Y,double* M,double* B,size_t n,size_t m)
{
	double d, *c;
	size_t i,j;
	size_t* idx = (size_t*)malloc(n*sizeof(size_t));
	if (!matrix_lu_double (M,idx,&d,n)) {
		free (idx);
		return 0;
	}
	c = (double*)malloc(n*sizeof(double));
	for (j=0;j<m;j++)
	{
		for (i=0;i<n;i++)
			c[i] = B[i+j*n];
		matrix_lub_double (c,M,idx,n);
		for (i=0;i<n;i++)
			Y[i+j*n] = c[i];
	}
	free (idx);
	free (c);
	return 1;
}

void rot_cols_double (double* M,size_t i,size_t k,size_t n,double c,double s)
{
	double temp;
	size_t j;
	for (j=0;j<n;j++)
  	{
    		temp = c *M[j+i*n] + s * M[j+k*n];
    		M[j+k*n] = -s * M[j+i*n] + c * M[j+k*n];
    		M[j+i*n] = temp;
  	}
}

int matrix_eigMax_double (double* l,double* v,double* M,size_t n)
{
	size_t i,j,k;
	double r;
	double* y = (double*)malloc(n*sizeof(double));
	*l = 1;
	for (i=0;i<n;i++)
		y[i] = 1.0/sqrt(n);;
	for (k=0;k<20;k++)
	{
		/* y = A'*v */
		for (i=0;i<n;i++)
		{
			r = 0;
			for (j=0;j<n;j++)
				r += v[j]*M[i+j*n];
			y[i] = r;
		}
		/* l = v'*y */
		*l = 0;
		for (i=0;i<n;i++)
			*l += v[i]*y[i];
		if (fabs(*l) < 1E-7) {
			free (y);
			return 0;
		}
		/* r = ||y - t*v||^2 */
		r = 0;
		for (i=0;i<n;i++)
			r += y[i] - *l*v[i];
		r = fabs(r/(*l));
		/*printf ("l=%f,r=%e\n",*l,r); */
		if (r < 1E-14)
			break;
		/* r = ||y||^2 */
		r = 0;
		for (i=0;i<n;i++)
			r += y[i]*y[i];
		if (fabs(r) < 1E-7) {
			free (y);
			return 0;
		}
		/* v = y/sqrt(r) */
		r = 1.0/sqrt(r);
		for (i=0;i<n;i++)
		{
			v[i] = y[i]*r;
		}
	}
	/* cherche la plus grande valeur de v; */
	/* si elle est negative, inverse v; */
	r = -1;
	for (i=0;i<n;i++)
		if (v[i]<r)
			r = v[i];
	if (r < 0)
	{
		for (i=0;i<n;i++)
			v[i] = -v[i];
	}
	free (v);
	return 1;
}

void matrix_eigSym_double (double* M,double* d,size_t n)
{
	double* e = (double*)malloc(n*sizeof(double));
	matrix_qtq_double (M,d,e,n);
	matrix_eigTrig_double (M,d,e,n);
	free (e);
}

void matrix_eigTrig_double (double* Q,double* d,double* e,size_t n)
{
	size_t i,i_min,i_max,split;
	double b_sqr, bk, ak1, bk1, ak2, bk2, z;
 	double c, c2, cs, s, s2, dd, mu;

	i_min = 0;
 	while (i_min < n)
 	{
    		i_max = n - 1;
    		for (i = i_min; i < n - 1; i++)
    		{
      			if (fabs(e[i]) < 1E-14)
      			{
				i_max = i;
				break;
      			}
    		}

    		if (i_max <= i_min)
    		{
      			i_min = i_max + 1;
      			continue;
    		}

		split = 0;
    		while (!split)
    		{
			dd = (d[i_max-1] - d[i_max]) / 2.0;
      			b_sqr = e[i_max-1] * e[i_max-1];
      			mu = d[i_max] - b_sqr / (dd + sgn(dd) * sqrt(dd * dd + b_sqr));

      			givens_double(d[i_min] - mu, e[i_min], &c, &s);
      			s = -s;
      			if (fabs(c) < SQRT2)
      			{
				c2 = c * c;
				s2 = 1.0 - c2;
      			}
      			else
      			{
				s2 = s * s;
				c2 = 1.0 - s2;
      			}
      			cs = c * s;
      			ak1 =	c2 * d[i_min] + s2 * d[i_min+1] - 2 * cs * e[i_min];
      			bk1 =	cs * (d[i_min] - d[i_min+1]) + (c2 - s2) * e[i_min];
      			ak2 =	s2 * d[i_min] + c2 * d[i_min+1] + 2 * cs * e[i_min];
      			bk2 = (i_min < i_max - 1) ? c * e[i_min+1] : 0.0;
      			z = (i_min < i_max - 1) ? -s * e[i_min+1] : 0.0;
      			d[i_min] = ak1;
      			d[i_min+1] = ak2;
      			e[i_min] = bk1;
      			if (i_min < i_max - 1)
				e[i_min+1] = bk2;
			rot_cols_double(Q, i_min, i_min+1, n, c, -s);

      			for (i = i_min + 1; i < i_max; i++)
      			{
				givens_double(e[i-1], z, &c, &s);
				s = -s;

				if (fabs(c) < SQRT2)
				{
	  				c2 = c * c;
	  				s2 = 1 - c2;
				}
				else
				{
		  			s2 = s * s;
		  			c2 = 1 - s2;
				}
				cs = c * s;
				bk = c * e[i-1] - s * z;
				ak1 = c2 * d[i] + s2 * d[i+1] - 2 * cs * e[i];
				bk1 = cs * (d[i] - d[i+1]) + (c2 - s2) * e[i];
				ak2 = s2 * d[i] + c2 * d[i+1] + 2 * cs * e[i];
				bk2 = (i + 1 < i_max) ? c * e[i+1] : 0.0;
				z = (i + 1 < i_max) ? -s * e[i+1] : 0.0;
				d[i] = ak1;
				d[i+1] = ak2;
				e[i] = bk1;
				/*printf ("%f;",e[i]); */
				if (i < i_max - 1)
	  				e[i+1] = bk2;
				if (i > i_min)
	  				e[i-1] = bk;
  				rot_cols_double(Q, i, i+1, n, c, -s);
      			}

      			for (i = i_min; i < i_max; i++)
      			{
				if (fabs(e[i]) < 1E-14)
				{
	  				e[i] = 0.0;
	  				split = 1;
				}
      			}
    		}
  	}
}

void matrix_eigTrig_values_double (double* d,double* e,size_t n)
{
	size_t i,i_min,i_max,split;
	double b_sqr, bk, ak1, bk1, ak2, bk2, z;
 	double c, c2, cs, s, s2, dd, mu;

	i_min = 0;
 	while (i_min < n)
 	{
    		i_max = n - 1;
    		for (i = i_min; i < n - 1; i++)
    		{
      			if (fabs(e[i]) < 1E-14)
      			{
				i_max = i;
				break;
      			}
    		}

    		if (i_max <= i_min)
    		{
      			i_min = i_max + 1;
      			continue;
    		}

		split = 0;
    		while (!split)
    		{
			dd = (d[i_max-1] - d[i_max]) / 2;
      			b_sqr = e[i_max-1] * e[i_max-1];
      			mu = d[i_max] - b_sqr / (dd + sgn(dd) * sqrt(dd * dd + b_sqr));

      			givens_double(d[i_min] - mu, e[i_min], &c, &s);
      			s = -s;
      			if (fabs(c) < SQRT2)
      			{
				c2 = c * c;
				s2 = 1 - c2;
      			}
      			else
      			{
				s2 = s * s;
				c2 = 1 - s2;
      			}
      			cs = c * s;
      			ak1 =	c2 * d[i_min] + s2 * d[i_min+1] - 2 * cs * e[i_min];
      			bk1 =	cs * (d[i_min] - d[i_min+1]) + (c2 - s2) * e[i_min];
      			ak2 =	s2 * d[i_min] + c2 * d[i_min+1] + 2 * cs * e[i_min];
      			bk2 = (i_min < i_max - 1) ? c * e[i_min+1] : 0.0;
      			z = (i_min < i_max - 1) ? -s * e[i_min+1] : 0.0;
      			d[i_min] = ak1;
      			d[i_min+1] = ak2;
      			e[i_min] = bk1;
      			if (i_min < i_max - 1)
				e[i_min+1] = bk2;

      			for (i = i_min + 1; i < i_max; i++)
      			{
				givens_double(e[i-1], z, &c, &s);
				s = -s;

				if (fabs(c) < SQRT2)
				{
	  				c2 = c * c;
	  				s2 = 1 - c2;
				}
				else
				{
		  			s2 = s * s;
		  			c2 = 1 - s2;
				}
				cs = c * s;
				bk = c * e[i-1] - s * z;
				ak1 = c2 * d[i] + s2 * d[i+1] - 2 * cs * e[i];
				bk1 = cs * (d[i] - d[i+1]) + (c2 - s2) * e[i];
				ak2 = s2 * d[i] + c2 * d[i+1] + 2 * cs * e[i];
				bk2 = (i + 1 < i_max) ? c * e[i+1] : 0.0;
				z = (i + 1 < i_max) ? -s * e[i+1] : 0.0;
				d[i] = ak1;
				d[i+1] = ak2;
				e[i] = bk1;
				/*printf ("%f;",e[i]); */
				if (i < i_max - 1)
	  				e[i+1] = bk2;
				if (i > i_min)
	  				e[i-1] = bk;
      			}

      			for (i = i_min; i < i_max; i++)
      			{
				if (fabs(e[i]) < 1E-14)
				{
	  				e[i] = 0.0;
	  				split = 1;
				}
      			}
    		}
  	}
}

void	matrix_sym_centering_double (double* M,size_t n,const double* buf)
{
    size_t i,j;
    double s = 0;
    double* b = NULL;
    if (buf == NULL)
	buf = b = (double*)malloc(n*sizeof(double));

    /* Somme des lignes/colonnes */
    for (j=0;j<n;j++) {
	double sum = 0;
	for (i=0;i<n;i++)
	    sum += M[i+j*n];
	s += sum;
	b[j] = sum / n;
    }
    s /= (n*n);

    /* Centre */
    for (j=0;j<n;j++) {
	double x = s - b[j];
	for (i=0;i<n;i++)
	    M[i+j*n] += x - b[i];
    }

    if (b)
	free(b);
}

size_t matrix_argmin_l2_double (const double* sample,const double* dict,size_t n,size_t m)
{
    size_t c,cMin = 0;
    float dMin = vector_l2p2_double(sample, dict, n);
    for (c = 1; c < m; c++) {
        float d = vector_l2p2_double(sample, dict+c*n, n);
        if (d < dMin) {
            dMin = d;
            cMin = c;
        }
    }
    return cMin;
}

int    matrix_save_double(const char* file_name,const double* M,size_t n,size_t m)
{
    size_t i,j;
    FILE* file = fopen(file_name,"wt");
    if (!file)
        return 0;
    for (j=0;j<m;j++) {
        for (i=0;i<n;i++) {
            if (i>0) fprintf(file,",");
            fprintf(file,"%f",M[i+j*n]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
    return 1;
}



void	matrix_CpAAt_double (double* C,const double* A,size_t n,size_t p)
{
	size_t i,j,k;
	for (k=0;k<p;k++)
	{
		for (j=0;j<n;j++)
		{
                        double r = A[j+k*n];
			for (i=0;i<n;i++)
				C[i+j*n] += r*A[i+k*n];
		}
	}
}



struct struct_sort_double
{
    size_t i;
    double val;
};

int compare_sort_double(const void * a, const void * b){
	double v = ((struct struct_sort_double *)b)->val - ((struct struct_sort_double *)a)->val;
    if (v > 0)
        return 1;
    if (v < 0)
        return -1;
    return 0;
}

void matrix_sortEig_double (double* A,double* d,size_t n)
{
    size_t i;
    struct struct_sort_double * sSort = (struct struct_sort_double *)malloc(n*sizeof(struct struct_sort_double));
    double * Atemp = (double *)malloc(n*n*sizeof(double));

    memcpy(Atemp, A, n*n*sizeof(double));

    for(i = 0 ; i < n ; i++)
    {
        sSort[i].i = i;
        sSort[i].val = d[i];
    }

    qsort(sSort, n, sizeof(struct struct_sort_double), compare_sort_double);

    for(i = 0 ; i < n ; i++)
    {
        d[i] = sSort[i].val;
        memcpy(&A[n*i], &Atemp[n*sSort[i].i], n*sizeof(double));
    }

    free(Atemp);
    free(sSort);
}

double	matrix_n2_double (const double* v,size_t n,size_t m)
{
	size_t i;
	double x;
        double s = 0;
	for (i=0;i<n*m;i++) {
            x = v[i];
            s += x*x;
        }
	return sqrt(s);
}
