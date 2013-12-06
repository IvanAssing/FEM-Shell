#include "matrix.h"


#include <lapacke.h>
#include <cblas.h>


#include <iomanip>
#include <cmath>

Matrix::Matrix()
{
}

Matrix::Matrix(int m_, int n_)
    :m(m_), n(n_)
{
    dim = m*n;
    data = new double[dim];
#pragma omp parallel for
    for(int i=0; i<dim; i++)
        data[i] = 0.0;
}

double& Matrix::operator()(int i, int j)
{
  return data[j*n + i];
}

const double& Matrix::operator()(int i, int j) const
{
  return data[j*n + i];
}

void Matrix::setUnit(int i)
{
    for(int k = 0; k<n ;k++)
        data[k*n + i] = 0.0;
    data[i*n + i] = 1.0;
}



void Matrix::operator = (Matrix &M)
{
    m = M.m;
    n = M.n;

    dim = m*n;
    data = new double[dim];
#pragma omp parallel for
    for(int i=0; i<dim; i++)
        data[i] = M.data[i];
}

void Matrix::operator_BtAB(Matrix &B) const
{
    Matrix C(m,m);

    for(int i=0; i<m; i++)
        for(int j=0; j<m; j++)
            for(int k=0; k<m; k++)
            C.data[j*m + i] += B.data[i*m + k]*data[j*m + k];

    for(int i=0; i<m; i++)
        for(int j=0; j<m; j++)
            for(int k=0; k<m; k++)
            data[j*m + i] += C.data[k*m + i]*B.data[j*m + k];
}




void Matrix::solve(Matrix &b, Matrix &x)
{

    int size = n;				/* dimension of matrix */


    int /*i,*/ j , c1, c2, Info;


    int *pivot = new int[size];



    c1=size;			/* and put all numbers we want to pass */
    c2=1;    			/* to the routine in variables */

    Info = LAPACKE_dgesv(LAPACK_COL_MAJOR, c1, c2, this->data, c1, pivot, b.data, c1);

    if(Info)
        std::cerr<<"Error in linear system solver: "<<Info;

    //    *  INFO    (output) INTEGER
    //    *          = 0:  successful exit
    //    *          < 0:  if INFO = -i, the i-th argument had an illegal value
    //    *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    //    *                has been completed, but the factor U is exactly
    //    *                singular, so the solution could not be computed.


    for (j=0; j<size; j++) x.data[j] = b.data[j];	/* print vector x */

    delete pivot;

}

std::ostream& operator<<(std::ostream& out, Matrix& obj)
{
    for(int i=0; i<obj.m; i++)
    {
        out<<"\n";
        for(int j=0; j<obj.n; j++)
            out<<std::setw(12)<<std::right<<obj.data[j*obj.n + i];
    }

    return out;
}


void Matrix::evalInverseMatrix(int n_, double **A)
{
    int info = 0;
    int*ipiv;
    double *_A;
    int dim = n_;

    _A = new double[dim*dim];
    ipiv = new int[dim];

    for(int i=0; i<n_ ;i++)
        for(int j=0; j<n_ ;j++)
            _A[i*dim+j] = A[i][j];

    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,  dim, dim, _A, dim, ipiv);

    if (info==0)
    {
        info = LAPACKE_dgetri(LAPACK_COL_MAJOR, dim, _A, dim, ipiv);
    }
    if(info)
        std::cerr<<"Error in matrix inverse";


    for(int i=0; i<n_ ;i++)
        for(int j=0; j<n_ ;j++)
            A[i][j] = _A[i*dim+j];

    delete [] ipiv;
    delete [] _A;
}

