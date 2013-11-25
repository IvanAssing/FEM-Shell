#include "matrix.h"

#include <f2c.h>
#include <clapack.h>

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

void Matrix::operator()(int i, int j, double value)
{
    data[j*n + i] = value;
}

void Matrix::add(int i, int j, double value)
{
    data[j*n + i] += value;
}

double Matrix::operator()(int i, int j)
{
    return data[j*n + i];
}

void Matrix::solve(Matrix &b, Matrix &x)
{
    //    void solveLinearSystem(unsigned n, double **A, double *b, double *x)
    //    {

    integer size = n;				/* dimension of matrix */


    integer /*i,*/ j , c1, c2, Info;


    //doublereal *AT = new doublereal[size*size];
    //doublereal *_b = new doublereal[size];
    integer *pivot = new integer[size];


    //        for (i=0; i<size; i++)
    //            _b[i] = b[i];

    //        for (i=0; i<size; i++)		/* to call a Fortran routine from C we */
    //        {				/* have to transform the matrix */
    //            for(j=0; j<size; j++) AT[j+size*i] = A[j][i];
    //        }

    c1=size;			/* and put all numbers we want to pass */
    c2=1;    			/* to the routine in variables */

    /* find solution using LAPACK routine SGESV, all the arguments have to */
    /* be pointers and you have to add an underscore to the routine name */

    dgesv_(&c1, &c2, this->data, &c1, pivot, b.data, &c1, &Info);

    //    *  INFO    (output) INTEGER
    //    *          = 0:  successful exit
    //    *          < 0:  if INFO = -i, the i-th argument had an illegal value
    //    *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
    //    *                has been completed, but the factor U is exactly
    //    *                singular, so the solution could not be computed.

    //#ifdef QT_GUI_LIB
    ////    int resp = QMessageBox::critical(this, QString("Erro na conectividade dos elementos"),
    ////                                     QString("O elemento %1 está conectado ao mesmo nó %2.\nVerifique a conectividade.").arg(elements[i]->index).arg(elements[i]->no1->index),
    ////                                     QMessageBox::Ok);


    //#endif


    /*
         parameters in the order as they appear in the function call
            order of matrix A, number of right hand sides (b), matrix A,
            leading dimension of A, array that records pivoting,
            result vector b on entry, x on exit, leading dimension of b
            return value */

    for (j=0; j<size; j++) x.data[j] = b.data[j];	/* print vector x */

    delete pivot;
    //        delete _b;
    //        delete AT;


    //    }


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
    integer info = 0;
    integer *ipiv;
    doublereal *_A;
    integer dim = n_;

    _A = new doublereal[dim*dim];
    ipiv = new integer[dim];

    for(int i=0; i<n_ ;i++)
        for(int j=0; j<n_ ;j++)
            _A[i*dim+j] = A[i][j];


    dgetrf_( &dim, &dim, _A, &dim, ipiv, &info );
    if (info==0)
    {
        doublereal workspace;
        integer tmp = -1;
        integer lwork;
        doublereal *work;
        dgetri_(&dim, _A, &dim, ipiv, &workspace, &tmp, &info);
        lwork = static_cast<int>(workspace);
        work = new doublereal[lwork];

        dgetri_( &dim, _A, &dim, ipiv, work, &lwork, &info );

        delete work;
    }
    if(info)
        std::cerr<<"Error in matrix inverse";



    for(int i=0; i<n_ ;i++)
        for(int j=0; j<n_ ;j++)
            A[i][j] = _A[i*dim+j];

    delete [] ipiv;
    delete [] _A;
}
