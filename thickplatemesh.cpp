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


    Matrix K(3*nNodes, 3*nNodes);

    Matrix f(3*nNodes, 1);
    Matrix x(3*nNodes, 1);

#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        dynamic_cast<ElementQN*>(elements[i])->getStiffnessMatrix(K, BftDBf, BctBc, L);



    #pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        if(nodes[i]->lockStatus[2])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 0, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 0, TP_NDOF*nodes[i]->index + 0) =  1.0;
        }
        if(nodes[i]->lockStatus[3])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 1, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 1, TP_NDOF*nodes[i]->index + 1) =  1.0;
        }
        if(nodes[i]->lockStatus[4])
        {
            for(int j=0; j<TP_NDOF*nNodes; j++)
                K(TP_NDOF*nodes[i]->index + 2, j) =  0.0;
            K(TP_NDOF*nodes[i]->index + 2, TP_NDOF*nodes[i]->index + 2) =  1.0;
        }
    }


    #pragma omp parallel for
    for(int i=0; i<nNodes; i++)
    {
        f(TP_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[2];
        f(TP_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[3];
        f(TP_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[4];
    }

    //    for(int i=0; i<ny; i++){
    //        for(int j=0; j<3*nNodes; j++)
    //        {
    //            K(3*nodes[i*nx]->index+0, j, 0.0);
    //            K(3*nodes[i*nx]->index+1, j, 0.0);
    //            K(3*nodes[i*nx]->index+2, j, 0.0);
    //        }
    //        K(3*nodes[i*nx]->index+0, 3*nodes[i*nx]->index+0, 1.0);
    //        K(3*nodes[i*nx]->index+1, 3*nodes[i*nx]->index+1, 1.0);
    //        K(3*nodes[i*nx]->index+2, 3*nodes[i*nx]->index+2, 1.0);
    //    }

    ////    for(int i=1; i<=ny; i++)
    ////        f(3*(i*nx-1), 0, -100.0/ny);

    //        f(3*nNodes-3, 0, -1200.);

    //        std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;

    //        std::cout<<"\n\n dim = "<<K.n<<std::flush;

    K.solve(f, x);

    std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;

    w = new double[nNodes];
    for(int i=0; i<nNodes; i++)
        w[i] = x(3*i, 0);

    //    //std::cout<<"\n\n TEMPO = "<<double(clock() - t_start)/CLOCKS_PER_SEC;

        for(int i=0; i<nNodes; i++)
            std::cout<<"\n"<<i<<"\t"<<x(3*i+0, 0)<<"\t"<<x(3*i+1, 0)<<"\t"<<x(3*i+2, 0);

    //plot(x);
}


void getMaxMin5(double *vector, int size, double &max, double &min)
{
    max = vector[0];
    min = vector[0];

    for(int i=1; i<size; i++){
        max = vector[i]>max ? vector[i] : max;
        min = vector[i]<min ? vector[i] : min;
    }
}

void ThickPlateMesh::draw(void)
{
    double T0, T1, T2, T3, T4, Tn, R, G, B;
    getMaxMin5(w, nNodes, T4, T0);

    T2 = (T0+T4)/2.;
    T1 = (T0+T2)/2.;
    T3 = (T2+T4)/2.;


    int k[4];

    for(int i=0; i<nElements; i++){

        int n = sqrt(elements[i]->np);


        k[0] = elements[i]->nodes[0]->index;
        k[1] = elements[i]->nodes[n-1]->index;
        k[2] = elements[i]->nodes[elements[i]->np-1]->index;
        k[3] = elements[i]->nodes[elements[i]->np-n]->index;

        glBegin(GL_QUADS);
        for(int p = 0; p<4; p++){
            Tn = w[k[p]];
            R = Tn<T2?  0. : (Tn>T3? 1. : (Tn-T2)/(T3-T2));
            B = Tn>T2?  0. : (Tn<T1? 1. : (T2-Tn)/(T2-T1));
            G = Tn<T1? (Tn-T0)/(T1-T0) : Tn>T3 ? (T4-Tn)/(T4-T3) : 1.;
            glColor4d(R,G,B,0.8);

            glVertex2d(nodes[k[p]]->x, nodes[k[p]]->y);
        }
        glEnd();

    }


}
