#include "elementssqn.h"

#include <QGLWidget>
#include "matrix.h"
#include "polynomial2d.h"
#include "integralgauss2d.h"
#include "rational2d.h"
#include <cmath>
#include <fstream>

#define S_NDOF 6

ElementSSQN::ElementSSQN(int npx_, int npy_, Node **nodes_, bool _selectiveIntegracion)
    :npx(npx_), npy(npy_), selectiveIntegration(_selectiveIntegracion)
{
    np = npx*npy;

    nodes = new Node*[np];

    for(int i=0; i<np; i++)
        nodes[i] = nodes_[i];
}

void ElementSSQN::draw(void)
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


void ElementSSQN::getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Polynomial2D **Bm, Lagrange *L)
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

    Matrix ke(S_NDOF*np, S_NDOF*np);

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    ke(S_NDOF*ii + 2 + i,  S_NDOF*ij + 2 + j) += IntegralGauss2D::intNP(npx, Bf[3*ii+i][3*ij+j]/J);

    int npi = selectiveIntegration? npx-1 : npx;

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<3; i++)
                for(int j=0; j<3; j++)
                    ke(S_NDOF*ii + 2 + i,  S_NDOF*ij + 2 + j) += IntegralGauss2D::intNP(npi, Bc[3*ii+i][3*ij+j]/J);

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<2; i++)
                for(int j=0; j<2; j++)
                    ke(S_NDOF*ii + i, S_NDOF*ij  + j) += IntegralGauss2D::intNP(npx, Bm[2*ii+i][2*ij+j]/J);


    // Matriz de rotacao
    double nx[3] = {nodes[1]->x - nodes[0]->x, nodes[1]->y - nodes[0]->y,  nodes[1]->z - nodes[0]->z};

    double len = sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2]);

    nx[0] /= len;
    nx[1] /= len;
    nx[2] /= len;

    double n_[3] = {nodes[npx]->x - nodes[0]->x, nodes[npx]->y - nodes[0]->y,  nodes[npx]->z - nodes[0]->z};

    double nz[3] = {nx[1]*n_[2] - nx[2]*n_[1], nx[2]*n_[0] - nx[0]*n_[2], nx[0]*n_[1] - nx[1]*n_[0]};

    len = sqrt(nz[0]*nz[0] + nz[1]*nz[1] + nz[2]*nz[2]);

    nz[0] /= len;
    nz[1] /= len;
    nz[2] /= len;

    double ny[3] = {nz[1]*nx[2] - nz[2]*nx[1], nz[2]*nx[0] - nz[0]*nx[2], nz[0]*nx[1] - nz[1]*nx[0]};

    Matrix R(S_NDOF*np, S_NDOF*np);

    for(int ii=0; ii<2*np; ii++)
        for(int ij=0; ij<2*np; ij++)
            for(int i=0; i<3; i++)
            {
                R(ii*3 + 0, ij*3 + i) = nx[i];
                R(ii*3 + 1, ij*3 + i) = ny[i];
                R(ii*3 + 2, ij*3 + i) = nz[i];
            }

    ke.operator_BtAB(R); // ke = Rt.ke.R

    for(int ii=0; ii<np; ii++)
        for(int ij=0; ij<np; ij++)
            for(int i=0; i<S_NDOF; i++)
                for(int j=0; j<S_NDOF; j++)
                    k(S_NDOF*nodes[ii]->index + i, S_NDOF*nodes[ij]->index + j) += ke(S_NDOF*ii + i, S_NDOF*ij  + j);

}



