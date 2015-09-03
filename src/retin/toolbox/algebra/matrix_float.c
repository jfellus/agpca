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

#include "matrix_float.h"

#define sgn(x) (((x) >= -1E-14) ? (1) : (-1))
#define SQRT2 (1.414213562373095048801688724209)


void	matrix_cpy_float (float* Y,size_t iy,size_t jy,const float* X,size_t ix,size_t jx,size_t n,size_t m,size_t nx,size_t ny)
{
    size_t i,j;
    for (j=0;j<m;j++) {
        for (i=0;i<n;i++) {
            Y[(i+iy)+(j+jy)*ny] = X[(i+ix)+(j+jx)*nx];
        }
    }
}

void	matrix_rot_float (float* Y,const float* X,size_t di,size_t dj,size_t n,size_t m)
{
    size_t i,j;
    for (j=0;j<m;j++) {
        size_t j0 = (j+dj)%m;
        for (i=0;i<n;i++) {
            size_t i0 = (i+di)%n;
            Y[i+j*n] = X[i0+j0*n];
        }
    }
}

void givens_float(float x, float y, float * c, float * s)
{
	float norm = sqrt(x * x + y * y);
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

void rot_cols_float (float* M,size_t i,size_t k,size_t n,float c,float s)
{
	float temp;
	size_t j;
	for (j=0;j<n;j++)
  	{
    		temp = c *M[j+i*n] + s * M[j+k*n];
    		M[j+k*n] = -s * M[j+i*n] + c * M[j+k*n];
    		M[j+i*n] = temp;
  	}
}


void	matrix_CpAB_float (float* pC,const float* pA,const float* pB,size_t n,size_t p,size_t m)
{
	float r;
	const float* p1;
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

void	matrix_CpAtB_float (float* pC,const float* pA,const float* pB,size_t n,size_t p,size_t m)
{
	size_t i/*,k*/,j;
	for (j=0;j<m;j++)
	{
		for (i=0;i<n;i++)
		{
			pC[i+j*n] += vector_ps_float(pA+i*p, pB+j*p, p);
		}
	}
}

void	matrix_CpABt_float (float* C,const float* A,const float* B,size_t n,size_t p,size_t m)
{
	float r;
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

#if defined(ALGEBRA_SSE2)

void	matrix_CpAAt_float (float* C,const float* A,size_t n,size_t p)
{
    size_t i,j,k;
    size_t q = n / 8;
    size_t r = n % 8;

    for (k=0;k<p;k++) {
        float* pC = C;
        for (j=0;j<n;j++) {
            __m128 w = _mm_load1_ps (A+j+k*n);
            const float* pA = A+k*n;
            if (ALGEBRA_IS_ALIGNED(pA) && ALGEBRA_IS_ALIGNED(pC)) {
                for (i=0;i<q;i++) {
                    __m128 i1 = _mm_load_ps(pA);
                    __m128 i2 = _mm_load_ps(pA+4);
                    __m128 o1 = _mm_load_ps(pC);
                    __m128 o2 = _mm_load_ps(pC+4);
                    _mm_store_ps(pC+0,_mm_add_ps(o1,_mm_mul_ps(i1,w)));
                    _mm_store_ps(pC+4,_mm_add_ps(o2,_mm_mul_ps(i2,w)));
                    pA += 8;
                    pC += 8;
                }
            }
            else {
                for (i=0;i<q;i++) {
                    __m128 i1 = _mm_loadu_ps(pA);
                    __m128 i2 = _mm_loadu_ps(pA+4);
                    __m128 o1 = _mm_loadu_ps(pC);
                    __m128 o2 = _mm_loadu_ps(pC+4);
                    _mm_storeu_ps(pC+0,_mm_add_ps(o1,_mm_mul_ps(i1,w)));
                    _mm_storeu_ps(pC+4,_mm_add_ps(o2,_mm_mul_ps(i2,w)));
                    pA += 8;
                    pC += 8;
                }
            }
            for (i=0;i<r;i++) {
                (*pC++) += A[j+k*n]*(*pA++);
            }
        }
    }
}

#else

void	matrix_CpAAt_float (float* C,const float* A,size_t n,size_t p)
{
	size_t i,j,k;
	for (k=0;k<p;k++)
	{
		for (j=0;j<n;j++)
		{
                        float r = A[j+k*n];
			for (i=0;i<n;i++)
				C[i+j*n] += r*A[i+k*n];
		}
	}
}

#endif

void	matrix_Cpaat_float (float* C,const float* a,size_t n)
{
	size_t j;
	for (j=0;j<n;j++)
		vector_addm_float(C+j*n, a[j], a, n);
}

void	matrix_CpABAt_float (float* C,const float* A,const float* B,size_t n,size_t m,float* r)
{
	float a;
	size_t i,k,t,j;
	for (j=0;j<n;j++)
	{
		memset (r,0,m*sizeof(float));
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

void	matrix_CpAtBA_float (float* C,const float* A,const float* B,size_t n,size_t m,float* r)
{
	float a;
	size_t i,k,t,j;
	for (j=0;j<n;j++)
	{
		memset (r,0,m*sizeof(float));
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


void	matrix_CpAdAt_float (float* C,const float* A,const float* d,size_t n,size_t m)
{
	float a;
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

void	matrix_t_float (float* Y,const float* X,size_t n,size_t m)
{
	size_t i,j;
	for (j=0;j<m;j++)
	for (i=0;i<n;i++)
		Y[j+i*m] = X[i+j*n];
}

void	matrix_sum_float (float* s,const float* X,size_t n,size_t m)
{
	float	r;
	size_t i,j;
	for (j=0;j<m;j++)
	{
		r = 0;
		for (i=0;i<n;i++)
			r += X[i+j*n];
		s[j] = r;
	}
}

void	matrix_sumt_float (float* s,const float* X,size_t n,size_t m)
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

void	matrix_meant_float (float* s,const float* X,size_t n,size_t m)
{
	size_t i,j;
        double* temp = (double*)malloc(n*sizeof(double));
        for (i=0;i<n;i++)
            temp[i] = 0;
	for (j=0;j<m;j++)
	{
            for (i=0;i<n;i++)
                    temp[i] += X[i+j*n];
	}
        for (i=0;i<n;i++)
            s[i] = temp[i] / ((double)m);
        free(temp);
}

#define matrix_float_gran 1024

double	matrix_l2p2_float (const float* v1,const float* v2,size_t n,size_t m)
{
	size_t i;
	double s = 0;
        size_t q = (n*m) / matrix_float_gran;
        size_t r = (n*m) % matrix_float_gran;
        for (i=0;i<q;i++)
            s += vector_l2p2_float(v1+i*matrix_float_gran,v2+i*matrix_float_gran,matrix_float_gran);
        s += vector_l2p2_float(v1+i*matrix_float_gran,v2+i*matrix_float_gran,r);
	return s;
}

double	matrix_ps_xxt_yyt_float (const float* x,const float* y,size_t n)
{
    double w = vector_ps_float(x,y,n);
    return w*w;
}

double	matrix_ps_XXt_YYt_float (const float* X,size_t p,const float* Y,size_t q,size_t n)
{
    size_t i,j;
    double sum = 0;
    for (j=0;j<q;j++) {
        for (i=0;i<p;i++) {
            double w = vector_ps_float(X+i*n,Y+j*n,n);
            sum += w*w;
        }
    }
    return sum;
}

double	matrix_n2_float (const float* v,size_t n,size_t m)
{
	size_t i;
	float x;
        double s = 0;
	for (i=0;i<n*m;i++) {
            x = v[i];
            s += x*x;
        }
	return sqrt(s);
}

double	matrix_l2_float (const float* pa,const float* pb,size_t n,size_t m)
{
    return sqrt(matrix_l2p2_float(pa,pb,n,m));
}



#define permuter(a,b) { r = (a); (a) = (b); (b) = r; }

size_t	matrix_ll_float (float* L,float* A,size_t n,size_t m)
{
	size_t i,j,k;
	float r;

	memcpy (L,A,n*m*sizeof(float));

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

void	matrix_inv_r_float (float* Y,float* L,size_t n)
{
    int i,j,k;
    float sum;
    memset(Y,0,n*n*sizeof(float));
    for (j=0;j<n;j++) {
        if (fabs(L[j+j*n]) < 1E-5)
            return;
        Y[j+j*n] = 1 / L[j+j*n];
        for (i=j-1;i>=0;i--) {
            for (sum=0,k=i+1;k<=j;k++) {
                sum -= L[i+k*n]*Y[k+j*n];
            }
            Y[i+j*n] = sum / L[i+i*n];
        }
    }
}

size_t	matrix_qr_float (float* pX,float* pR,size_t n,size_t m)
{
	size_t i,j,k,rang;
	float r,x;
	memset (pR,0,m*m*sizeof(float));
	for (j=0;j<m;j++)
	{
		r = 0;
		for (i=0;i<n;i++)
		{
			x = pX[i+j*n];
			r += x*x;
		}
		if (r < 1E-5)
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
        rang = j;
	for (;j<m;j++)
	{
		for (i=0;i<n;i++)
			pX[i+j*n] = 0;
	}
        return rang;
}

size_t	matrix_inv_cov_float (float* invC,float* X,size_t n,size_t m)
{
    float* R;
    size_t rang;
    
    if (n < m)
        return 0;

    R = (float*)malloc(m*m*sizeof(float));
    rang = matrix_qr_float(X,R,n,m);
    matrix_inv_r_float (X,R,m);
    memset(invC,0,m*m*sizeof(float));
    matrix_CpABt_float (invC,X,X,m,m,m);

    free(R);
    return rang;
}

void	matrix_lle_float(float* y,const float* x,const float* X,const float* invC,size_t n,size_t m)
{
    size_t i,k;

   /* double alpha = 1;
    double beta = 0;
    for (i=0;i<m;i++) {
	for (k=0;k<m;k++) {
	    double t = invC[k+i*m];
	    beta += t;
	    alpha -= t*vector_ps_float(x,X+k*n,n);
	}
    }
    double lambda = alpha / beta;
    printf ("lambda = %f\n",lambda);*/
    for (i=0;i<m;i++) {
	double sum = 0;
	for (k=0;k<m;k++) {
	    double t = invC[k+i*m];
	    sum += t*vector_ps_float(x,X+k*n,n); 
	}
	y[i] = sum;
    }
}

void	matrix_qtq_float (float* A,float* d,float* e,size_t n)
{
	size_t i,j,k;
	float rho,beta;
	float* r = (float*)malloc(n*sizeof(float));
	float* s = (float*)malloc(n*sizeof(float));
	memset (r,0,n*sizeof(float));
	memset (s,0,n*sizeof(float));
	for (k=0;k<n-2;k++)
	{
		d[k] = A[k+k*n];
		beta = 0.0;
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

void	matrix_qtq_t_float (float* A,float* d,float* e,size_t n)
{
	size_t i,j,k;
	float rho,beta;
	float* r = (float*)malloc(n*sizeof(float));
	float* s = (float*)malloc(n*sizeof(float));
	memset (r,0,n*sizeof(float));
	memset (s,0,n*sizeof(float));
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

int	matrix_lu_float (float* M,size_t* idx,float* d,size_t n)
{
	size_t i,imax,j,k;
	float big,dum,sum,temp;
	float* v = (float*)malloc(n*sizeof(float));
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

void matrix_lub_float (float* b,float* M,size_t *idx,size_t n)
{
	int i,ii=-1,ip,j;
	float sum;
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

int matrix_inv_float (float* Y,float* M,size_t n)
{
	float d, *c;
	size_t i,j;
	size_t* idx = (size_t*)malloc(n*sizeof(size_t));
	if (!matrix_lu_float (M,idx,&d,n)) {
		free (idx);
		return 0;
	}
	c = (float*)malloc(n*sizeof(float));
	for (j=0;j<n;j++)
	{
		for (i=0;i<n;i++)
			c[i] = 0.0;
		c[j] = 1.0;
		matrix_lub_float (c,M,idx,n);
		for (i=0;i<n;i++)
			Y[i+j*n] = c[i];
	}
	free (idx);
	free (c);
	return 1;
}

int matrix_invm_float (float* Y,float* M,float* B,size_t n,size_t m)
{
	float d, *c;
	size_t i,j;
	size_t* idx = (size_t*)malloc(n*sizeof(size_t));
	if (!matrix_lu_float (M,idx,&d,n)) {
		free (idx);
		return 0;
	}
	c = (float*)malloc(n*sizeof(float));
	for (j=0;j<m;j++)
	{
		for (i=0;i<n;i++)
			c[i] = B[i+j*n];
		matrix_lub_float (c,M,idx,n);
		for (i=0;i<n;i++)
			Y[i+j*n] = c[i];
	}
	free (idx);
	free (c);
	return 1;
}

int matrix_eigMax_float (float* l,float* v,float* M,size_t n)
{
	size_t i,j,k;
	float r;
	float* y = (float*)malloc(n*sizeof(float));
	*l = 1;
	for (i=0;i<n;i++)
		v[i] = 1.0/sqrt(n);;
	for (k=0;k<50;k++)
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
	free (y);
	return 1;
}

int matrix_eigMaxSym_float (float* l,float* v,float* M,size_t n){
    size_t i,k;
    float r;
    float* y = (float*)malloc(n*sizeof(float));
    *l = 1;
    for (i=0;i<n;i++)
            v[i] = 1.0/sqrt(n);
    for (k=0;k<500;k++)
    {
        /* y = A'*v */
        for (i=0;i<n;i++)
        {
            y[i] = vector_ps_float (v, &M[i*n], n);
        }
        /* l = v'*y */
        *l = vector_ps_float (v, y, n);
        if (fabs(*l) < 1E-7) {
            free (y);
            return 0;
        }
        /* r = ||y - t*v||^2 */
        vector_linear_float (v,y,-(*l),v,n);
        r = vector_ps_float (v, v, n);
        r = fabs(r/(*l));
        /*printf ("l=%f,r=%e\n",*l,r); */
        if (r < 1E-14){
            vector_cpym_float (v, 1.0/vector_n2_float (y, n), y, n);
            break;
        }
        /* r = ||y||^2 */
        r = vector_n2p2_float (y, n);
        if (r < 1E-7) {
            free (y);
            return 0;
        }
        /* v = y/sqrt(r) */
        r = 1.0/sqrt(r);
        vector_cpym_float (v, r, y, n);
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
    free (y);
    return 1;
}

void matrix_eigSym_float (float* M,float* d,size_t n)
{
	float* e = (float*)malloc(n*sizeof(float));
	matrix_qtq_float (M,d,e,n);
	matrix_eigTrig_float (M,d,e,n);
	free (e);
}

struct struct_sort_float
{
    size_t i;
    float val;
};

int compare_sort_float(const void * a, const void * b){
    float v = ((struct struct_sort_float *)b)->val - ((struct struct_sort_float *)a)->val;
    if (v > 0)
        return 1;
    if (v < 0)
        return -1;
    return 0;
}

void matrix_sortEig_float (float* A,float* d,size_t n)
{
    size_t i;
    struct struct_sort_float * sSort = (struct struct_sort_float *)malloc(n*sizeof(struct struct_sort_float));
    float * Atemp = (float *)malloc(n*n*sizeof(float));
    
    memcpy(Atemp, A, n*n*sizeof(float));
    
    for(i = 0 ; i < n ; i++)
    {
        sSort[i].i = i;
        sSort[i].val = d[i];
    }
    
    qsort(sSort, n, sizeof(struct struct_sort_float), compare_sort_float);
    
    for(i = 0 ; i < n ; i++)
    {
        d[i] = sSort[i].val;
        memcpy(&A[n*i], &Atemp[n*sSort[i].i], n*sizeof(float));
    }
    
    free(Atemp);
    free(sSort);
}

void matrix_eigTrig_float (float* Q,float* d,float* e,size_t n)
{
	size_t i,i_min,i_max,split;
	float b_sqr, bk, ak1, bk1, ak2, bk2, z;
 	float c, c2, cs, s, s2, dd, mu;

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

      			givens_float(d[i_min] - mu, e[i_min], &c, &s);
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
			rot_cols_float(Q, i_min, i_min+1, n, c, -s);

      			for (i = i_min + 1; i < i_max; i++)
      			{
				givens_float(e[i-1], z, &c, &s);
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
				/* printf ("%f;",e[i]); */
				if (i < i_max - 1)
	  				e[i+1] = bk2;
				if (i > i_min)
	  				e[i-1] = bk;
  				rot_cols_float(Q, i, i+1, n, c, -s);
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

void matrix_eigTrig_values_float (float* d,float* e,size_t n)
{
	size_t i,i_min,i_max,split;
	float b_sqr, bk, ak1, bk1, ak2, bk2, z;
 	float c, c2, cs, s, s2, dd, mu;

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

      			givens_float(d[i_min] - mu, e[i_min], &c, &s);
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
				givens_float(e[i-1], z, &c, &s);
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

void	matrix_sym_centering_float (float* M,size_t n,const float* buf)
{
    size_t i,j;
    float s = 0;
    float* b = NULL;
    if (buf == NULL)
	buf = b = (float*)malloc(n*sizeof(float));    
    
    /* Somme des lignes/colonnes */
    for (j=0;j<n;j++) {
	float sum = 0;
	for (i=0;i<n;i++)
	    sum += M[i+j*n];
	s += sum;
	b[j] = sum / n;
    }
    s /= (n*n);

    /* Centre */
    for (j=0;j<n;j++) {
	float x = s - b[j];
	for (i=0;i<n;i++)
	    M[i+j*n] += x - b[i];
    }

    if (b)
	free(b);
}

double matrix_cps_float(const float* K,const float* L,size_t n)
{
    size_t i,j;
    double a = 0;
    double sumL = 0;
    double sumK = 0;
    for (j=0;j<n;j++) {
        float s = 0;
	float colL = 0;
	float colK = 0;
	for (i=0;i<n;i++) {
            float x1 = K[i+j*n];
            float x2 = L[i+j*n];
	    colK += x1;
	    colL += x2;
	    s += x1 * x2;
	}
	sumL += colL;
	sumK += colK;
	a += ((double)s) - 2*((double)colL)*((double)colK) / n;
    }
    return a + sumL*sumK / (n*n);
}

double matrix_cn2_float (const float* K,size_t n)
{
    size_t i,j;
    double a = 0;
    double sum = 0;
    for (j=0;j<n;j++) {
        float s = 0;
	float col = 0;
	for (i=0;i<n;i++) {
	    float x = K[i+j*n];
	    col += x;
	    s += x*x;
	}
	sum += col;
	a += ((double)s) - 2*((double)col)*((double)col) / n;
    }
    return a + sum*sum / (n*n);
}

double	matrix_cps_xxt_yyt_float (const float* x,const float* y,size_t n)
{
    size_t i;
    float ps = 0;
    float s1 = 0, s2 = 0;
    double w;
    for (i=0;i<n;i++) {
        float x1 = x[i];
        float x2 = y[i];
        ps += x1*x2;
        s1 += x1;
        s2 += x2;
    }
    w = ps - ((double)s1)*((double)s2)/n;
    return w*w;
}

/* en matlab:
  function [r] = cps_XXt_yyt (X,y)
    A=X'*(y - sum(y)/n);
    r=sum(A.*A);
  endfunction
 */
double	matrix_cps_XXt_yyt_float (const float* X,size_t p,const float* y,size_t n)
{
    size_t i,j;
    float ps, s1;
    double m2 = 0, sum = 0, w;
    for (i=0;i<n;i++)
        m2 += y[i];
    m2 /= n;
    for (j=0;j<p;j++) {
        ps = 0;
        s1 = 0;
        for (i=0;i<n;i++) {
            float x1 = X[i+j*n];
            ps += x1*y[i];
            s1 += x1;
        }
        w = ps - s1*m2;
        sum += w*w;
    }
    return sum;
}

void matrix_diag_XXt_yyt_float (const float* X,size_t p,const float* y,size_t n,float *res)
{
    size_t i,j;
    double m1 = 0, m2, sum;
    for (i=0;i<n;i++){
        m1 += y[i];
        res[i]=0;
    }
    m1 /= n;
    for (j=0;j<p;j++) {
        m2=0;
        for (i=0;i<n;i++) 
            m2 += X[i+j*n];        
        m2 /=n;
        sum=0;
        for (i=0;i<n;i++) 
            sum+=(y[i]-m1)*(X[i+j*n]-m2);        
        for (i=0;i<n;i++)
            res[i]+=X[i+j*n]*sum;
    }

    for (i=0;i<n;i++) {
        res[i]*=(y[i]-m1);
    }
}

/* en matlab:

function [r] = cps_XXt_YYt (X,Y)
  A=X'*Y - sum(X,1)'*sum(Y,1)/n;
  r=sum(A.*A);
endfunction

 */

double	matrix_cps_XXt_YYt_float (const float* X,size_t p,const float* Y,size_t q,size_t n)
{
    size_t i,r,s;
    float ps, s1;
    double sum = 0, w;
    for (s=0;s<q;s++) {
        double m2 = 0;
        for (i=0;i<n;i++)
            m2 += Y[i+s*n];
        m2 /= n;
        for (r=0;r<p;r++) {
            ps = 0;
            s1 = 0;
            for (i=0;i<n;i++) {
                float x1 = X[i+r*n];
                float x2 = Y[i+s*n];
                ps += x1*x2;
                s1 += x1;
            }
            w = ps - s1*m2;
            sum += w*w;
        }
    }
    return sum;
}

double	matrix_cps_XtX_YtY_float (const float* X,size_t p,const float* Y,size_t q,size_t n)
{

    size_t i,r,s;
    double sum = 0;
    double sum1 = 0;
    double sum2 = 0;

    for (r=0;r<p;r++) {
        double c1 = 0;
        for (i=0;i<n;i++)
            c1 += X[r+i*p];
        for (s=0;s<q;s++) {
            double c2 = 0;
            double ps = 0;
            for (i=0;i<n;i++) {
                double x = X[r+i*p];
                double y = Y[s+i*q];
                ps += x*y;
                c2 += y;
            }
            sum += ps*ps - 2*c1*c2*ps / n;
        }
        sum1 += c1*c1;
       
    }
    for (s=0;s<q;s++) {
        double c2 = 0;
        for (i=0;i<n;i++)
            c2 += Y[s+i*q];
         sum2 += c2*c2;
    }
    return sum + sum1*sum2 / (n*n);


/*    size_t i,j,r,s;
    double sum = 0;
    double sum1 = 0;
    double sum2 = 0;
    for (j=0;j<n;j++) {
	double col1 = 0;
	double col2 = 0;
	for (i=0;i<n;i++) {
            double x = 0;
            for (r=0;r<p;r++)
                x += X[r+i*p]*X[r+j*p];

            double y = 0;
            for (s=0;s<q;s++)
                y += Y[s+i*q]*Y[s+j*q];

	    sum += x*y;
	    col1 += x;
	    col2 += y;
	}
	sum1 += col1;
	sum2 += col2;
	sum = sum - 2*col1*col2 / n;
    }
    return sum + sum1*sum2 / (n*n);*/
}

size_t matrix_argmin_l2_float (const float* sample,const float* dict,size_t n,size_t m)
{
    size_t c,cMin = 0;
    float dMin = vector_l2p2_float(sample, dict, n);
    for (c = 1; c < m; c++) {
        float d = vector_l2p2_float(sample, dict+c*n, n);
        if (d < dMin) {
            dMin = d;
            cMin = c;
        }
    }
    return cMin;
}

void matrix_selft_float(float *data,size_t n,size_t m)
{
    float temp,swap;
    size_t abort;
    size_t i, j, k, k_start, k_new;
    size_t length = n*m;
 
    for(k_start=1; k_start < length; k_start++)
    {
        temp    = data[k_start];
        abort   = 0;
        k_new = k = k_start;
        do {
            if( k_new < k_start )
            {
                abort = 1;
                break;
            }
 
            k       = k_new;
            i       = k%n;
            j       = k/n;
            k_new   = i*m + j;
        }while(k_new != k_start);
 
 
        // if the current value is not the minimum of the cycle, then don't
        // perform the cycle
        if(abort)
            continue;
 
 
        // otherwise, perform the cycle
        k_new = k = k_start;
        do
        {
            swap = temp;
            temp = data[k_new];
            data[k_new] = swap;
            
            k       = k_new;
            i       = k%n;
            j       = k/n;
            k_new   = i*m + j;
        }while(k_new != k_start);
 
        swap = temp;
        temp = data[k_new];
        data[k_new] = swap;
    }
}

