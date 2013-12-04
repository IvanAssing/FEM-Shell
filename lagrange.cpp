#include "lagrange.h"
#include "matrix.h"

Lagrange::Lagrange(int n1_, int n2_)
    :n1(n1_), n2(n2_)
{
    int nn_ = (n1+1)*(n2+1);
    ne = nn_;
    double **a = new double*[nn_];
     for(int i=0; i<nn_; i++)
         a[i] = new double[nn_];

     int npx = n1+1, npy = n2+1;
     double dx = 2./static_cast<double>(npx-1);
     double dy = 2./static_cast<double>(npy-1);
     double x,y,xj,yi;

     for(int ny=0; ny<npy; ny++)
     {
         y = -1.+dy*ny;
         for(int nx=0; nx<npx; nx++)
         {
             x = -1.+dx*nx;
             yi = 1.;

             a[ny*npy+nx][0] = 1.;
             xj = 1.;
             for(int j=1; j<npx; j++)
             {
                 xj *= x;
                 a[ny*npy+nx][j] = xj;
             }

             for(int i=1; i<npy; i++)
             {
                 yi *= y;
                 xj = 1.;
                 a[ny*npy+nx][i*npy] = yi;
                 for(int j=1; j<npx; j++)
                 {
                     xj *= x;
                     a[ny*npy+nx][i*npy+j] = yi*xj;
                 }
             }
         }
     }

     Matrix::evalInverseMatrix(nn_, a);

     double *an = new double[nn_];
     N = new Polynomial2D[nn_];
     D1 = new Polynomial2D[nn_];
     D2 = new Polynomial2D[nn_];

     for(int i=0; i<nn_; i++){
         for(int j=0; j<nn_; j++)
             an[j] = a[j][i];
         N[i] = Polynomial2D(n1, n2, an);
         D1[i] = N[i].differential1();
         D2[i] = N[i].differential2();
     }

     delete [] an;

     for(int i=0; i<nn_; i++)
         delete [] a[i];

     delete [] a;
}

Lagrange::~Lagrange()
{
    delete [] D2;
    delete [] D1;
    delete [] N;

}
