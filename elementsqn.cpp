#include "elementsqn.h"

#include <QGLWidget>
#include "matrix.h"
#include "polynomial2d.h"
#include "integralgauss2d.h"

#include <cmath>



ElementSQN::ElementSQN(int npx_, int npy_, Node **nodes_, bool _selectiveIntegracion)
    :npx(npx_), npy(npy_), selectiveIntegration(_selectiveIntegracion)
{
    np = npx*npy;

    nodes = new Node*[np];

    for(int i=0; i<np; i++)
        nodes[i] = nodes_[i];
}

void ElementSQN::draw(void)
{
    glColor4d(0.0, 1.0, 0.0, 0.8);
    glLineWidth(2.5f);

    int n = npx;

    glBegin(GL_LINE_LOOP);
    {
        glVertex3d(nodes[0]->x, nodes[0]->y, nodes[0]->z);
        glVertex3d(nodes[n-1]->x, nodes[n-1]->y, nodes[n-1]->z);
        glVertex3d(nodes[np-1]->x, nodes[np-1]->y, nodes[np-1]->z);
        glVertex3d(nodes[np-n]->x, nodes[np-n]->y, nodes[np-n]->z);
    }
    glEnd();

}


void ElementSQN::getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Polynomial2D **Bm, Lagrange *L)
{

    Polynomial2D J, dxd1, dxd2, dyd1, dyd2;

    for(int i=0; i<np; i++)
    {
        dxd1 = dxd1 + L->D1[i]*nodes[i]->x;
        dxd2 = dxd2 + L->D2[i]*nodes[i]->x;
        dyd1 = dyd1 + L->D1[i]*nodes[i]->y;
        dyd2 = dyd2 + L->D2[i]*nodes[i]->y;
    }

    J = dxd1*dyd2 - dxd2*dyd1;


    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    k(5*nodes[ii]->index + 2 + i, 5*nodes[ij]->index + 2 + j) += IntegralGauss2D::intNP(npx, J * Bf[3*ii+i][3*ij+j]);
    //k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::int10P(J * Bf[3*ii+i][3*ij+j]);

    int npi = selectiveIntegration? npx-1 : npx;

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    k(5*nodes[ii]->index + 2 + i , 5*nodes[ij]->index + 2 + j) += IntegralGauss2D::intNP(npi, J * Bc[3*ii+i][3*ij+j]);
    //k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::int10P(J * Bc[3*ii+i][3*ij+j]);

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<2; i++)
                for(int j=0; j<2; j++)
                    k(5*nodes[ii]->index + i, 5*nodes[ij]->index  + j) += IntegralGauss2D::intNP(npx, J * Bm[2*ii+i][2*ij+j]);


}


