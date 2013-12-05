#include "thickplatemesh.h"

#include <QGLWidget>
#include <cmath>

#define TP_NDOF 3

ThickPlateMesh::ThickPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementQN **_elements, int _npx, int _npy, Matrix _D, double _GKt)
    :nNodes(_nNodes), nodes(_nodes), nElements(_nElements), elements(_elements), npx(_npx+1), npy(_npy+1), D(_D), GKt(_GKt)
{
    L = new Lagrange(npx-1, npy-1);

    npt = TP_NDOF*npx*npy;

    Polynomial2D Bf[3][npt];
    Polynomial2D Bc[2][npt];

    for(int i=0; i<npx*npy; i++)
    {
        Bf[0][3*i+2] = -1.0*L->D1[i];
        Bf[1][3*i+1] = L->D2[i];
        Bf[2][3*i+1] = L->D1[i];
        Bf[2][3*i+2] = -1.0*L->D2[i];

        Bc[0][3*i] = L->D1[i];
        Bc[0][3*i+2] = L->N[i];
        Bc[1][3*i] = L->D2[i];
        Bc[1][3*i+1] = -1.0*L->N[i];
    }

    BftDBf = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BftDBf[i] = new Polynomial2D[npt];

    Polynomial2D DB[3][npt];

    for(int i=0; i<3; i++)
        for(int j=0; j<npt; j++)
            DB[i][j] = Bf[0][j]*D(i,0) + Bf[1][j]*D(i,1) + Bf[2][j]*D(i,2);

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BftDBf[i][j] = Bf[0][i]*DB[0][j] + Bf[1][i]*DB[1][j] + Bf[2][i]*DB[2][j];

    BctBc = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BctBc[i] = new Polynomial2D[npt];

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BctBc[i][j] = (Bc[0][i]*Bc[0][j] + Bc[1][i]*Bc[1][j])*GKt;

}



void ThickPlateMesh::solve(void)
{

    int sys_dim = TP_NDOF*nNodes;

    Matrix K(sys_dim , sys_dim );


#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->getStiffnessMatrix(K, BftDBf, BctBc, L);


#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        if(nodes[i]->lockStatus[2])
            K.setUnit(TP_NDOF*nodes[i]->index + 0);
        if(nodes[i]->lockStatus[3])
            K.setUnit(TP_NDOF*nodes[i]->index + 1);
        if(nodes[i]->lockStatus[4])
            K.setUnit(TP_NDOF*nodes[i]->index + 2);
    }


    Matrix f(sys_dim);


#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        f(TP_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[2];
        f(TP_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[3];
        f(TP_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[4];
    }

    Matrix x(sys_dim);


    K.solve(f, x);

    results = new double*[2*TP_NDOF];
    for(int i=0; i<2*TP_NDOF; i++)
        results[i] = new double[nNodes];

#pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        results[0][i] = x(TP_NDOF*i + 0);
        results[1][i] = x(TP_NDOF*i + 1);
        results[2][i] = x(TP_NDOF*i + 2);
    }


}


void ThickPlateMesh::draw(vout option)
{
    double *x;
    switch (option) {
    case W:
        x = results[0];
        break;
    case RX:
        x = results[1];
        break;
    case RY:
        x = results[2];
        break;
    default:
        break;
    }

    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin(x, nNodes, T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;


    int k[4];

    for(int i=0; i<nElements; i++){

        glBegin(GL_QUADS);
        for(int ii = 0; ii<npy-1; ii++)
            for(int jj = 0; jj<npx-1; jj++)
            {
                k[0] = elements[i]->nodes[ii*npx+jj]->index;
                k[1] = elements[i]->nodes[ii*npx+jj+1]->index;
                k[2] = elements[i]->nodes[(ii+1)*npx+jj+1]->index;
                k[3] = elements[i]->nodes[(ii+1)*npx+jj]->index;

                for(int p = 0; p<4; p++){
                    Tn = x[k[p]];
                    R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
                    B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
                    G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
                    glColor4d(R,G,B,0.8);

                    glVertex2d(nodes[k[p]]->x, nodes[k[p]]->y);
                }
            }
        glEnd();

    }

    for(int i=0; i<nElements; i++)
        elements[i]->draw();

    for(int i=0; i<nNodes; i++)
        nodes[i]->draw_lock();

    for(int i=0; i<nNodes; i++)
        nodes[i]->draw_load();

    for(int i=0; i<nNodes; i++)
        nodes[i]->draw();

}
