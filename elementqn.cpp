#include "elementqn.h"
#include <QGLWidget>
#include "matrix.h"
#include "polynomial2d.h"
#include "integralgauss2d.h"

#include <cmath>


ElementQN::ElementQN(int npx_, int npy_, Node **nodes_, bool _selectiveIntegracion)
    :npx(npx_), npy(npy_), selectiveIntegration(_selectiveIntegracion)
{
    np = npx*npy;
    nodes = new Node*[np];

    for(int i=0; i<np; i++)
        nodes[i] = nodes_[i];
}

void ElementQN::draw(void)
{

    glColor4d(0.0, 1.0, 0.0, 0.8);
    glLineWidth(2.5f);

    glBegin(GL_LINE_LOOP);
    {
        glVertex3d(nodes[0]->x, nodes[0]->y, nodes[0]->z);
        glVertex3d(nodes[npx-1]->x, nodes[npx-1]->y, nodes[npx-1]->z);
        glVertex3d(nodes[np-1]->x, nodes[np-1]->y, nodes[np-1]->z);
        glVertex3d(nodes[np-npx]->x, nodes[np-npx]->y, nodes[np-npx]->z);
    }
    glEnd();
}


void ElementQN::getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Lagrange *L)
{

    Polynomial2D dxd1, dxd2, dyd1, dyd2;

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
                    k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::intNP(npx, J * Bf[3*ii+i][3*ij+j]);
    //k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::int10P(J * Bf[3*ii+i][3*ij+j]);

    int npi = selectiveIntegration? npx-1 : npx;

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::intNP(npi, J * Bc[3*ii+i][3*ij+j]);
    //k(3*nodes[ii]->index + i, 3*nodes[ij]->index  + j) += IntegralGauss2D::int10P(J * Bc[3*ii+i][3*ij+j]);


}

void ElementQN::evalResults(Matrix &M, Matrix &Q, Matrix &U, Polynomial2D **Bf, Polynomial2D **Bc)
{
    Polynomial2D DBf[3][3*np];
    Polynomial2D DBc[2][3*np];

    for(int i=0; i<3; i++)
        for(int j=0; j<3*np; j++)
            DBf[i][j] = J*Bf[i][j];


    for(int i=0; i<2; i++)
        for(int j=0; j<3*np; j++)
            DBc[i][j] = J*Bc[i][j];

    Matrix Me(3*np);

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    Me(3*ii+i) += DBf[i][3*ij + j](nodes[ii]->x, nodes[ii]->y)*U(3*nodes[ij]->index + j);

    Matrix Qe(2*np);

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<2; i++)
                for(int j=0; j<3; j++)
                    Qe(2*ii+i) += DBc[i][3*ij + j](nodes[ii]->x, nodes[ii]->y)*U(3*nodes[ij]->index + j);


    for(int ii=0; ii<np; ii++)
        for(int i=0; i<3; i++)
            M(3*nodes[ii]->index + i) = Me(3*ii+i);

    for(int ii=0; ii<np; ii++)
        for(int i=0; i<2; i++)
            Q(2*nodes[ii]->index + i) = Qe(2*ii+i);

}


