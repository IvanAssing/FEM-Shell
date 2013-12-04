#include "flatshellmesh.h"

#define TP_NDOF 3
#define PS_NDOF 2
#define FS_NDOF 5

FlatShellMesh::FlatShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSQN **_elements, int _npx, int _npy, Matrix _Df, double _GKt, Matrix _Dm)
    :nNodes(_nNodes), nodes(_nodes), nElements(_nElements), elements(_elements), npx(_npx+1), npy(_npy+1), Df(_Df), GKt(_GKt), Dm(_Dm)
{

    L = new Lagrange(npx-1, npy-1);

        npt = TP_NDOF*npx*npy;
        int npps = PS_NDOF*npx*npy;

    Polynomial2D Bf[3][npt];
    Polynomial2D Bc[2][npt];
    Polynomial2D Bm[3][npps];

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

        Bm[0][2*i] = L->D1[i];
        Bm[1][2*i+1] = L->D2[i];
        Bm[2][2*i] = L->D2[i];
        Bm[2][2*i+1] = L->D1[i];
    }


    BftDBf = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BftDBf[i] = new Polynomial2D[npt];

    Polynomial2D DB[3][npt];

    for(int i=0; i<3; i++)
        for(int j=0; j<npt; j++)
            DB[i][j] = Bf[0][j]*Df(i,0) + Bf[1][j]*Df(i,1) + Bf[2][j]*Df(i,2);

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BftDBf[i][j] = Bf[0][i]*DB[0][j] + Bf[1][i]*DB[1][j] + Bf[2][i]*DB[2][j];


    BctBc = new Polynomial2D*[npt];
    for(int i=0; i<npt; i++)
        BctBc[i] = new Polynomial2D[npt];

    for(int i=0; i<npt; i++)
        for(int j=0; j<npt; j++)
            BctBc[i][j] = (Bc[0][i]*Bc[0][j] + Bc[1][i]*Bc[1][j])*GKt;



    // *********************


    BmtDBm = new Polynomial2D*[npps];
    for(int i=0; i<npps; i++)
        BmtDBm [i] = new Polynomial2D[npps];

    Polynomial2D DmB[3][npps];

    for(int i=0; i<3; i++)
        for(int j=0; j<npps; j++)
            DmB[i][j] = Bm[0][j]*Dm(i,0) + Bm[1][j]*Dm(i,1) + Bm[2][j]*Dm(i,2);

    for(int i=0; i<npps; i++)
        for(int j=0; j<npps; j++)
            BmtDBm[i][j] = Bm[0][i]*DmB[0][j] + Bm[1][i]*DmB[1][j] + Bm[2][i]*DmB[2][j];


    Matrix K(5*nNodes, 5*nNodes);

    Matrix f(5*nNodes);
    Matrix x(5*nNodes);

#pragma omp parallel for
    for(int i=0; i<nElements; i++)
        elements[i]->getStiffnessMatrix(K, BftDBf, BctBc, BmtDBm, L);



#pragma omp parallel for
for(int i=0; i<nNodes; i++)
{
    if(nodes[i]->lockStatus[0])
    {
        for(int j=0; j<FS_NDOF*nNodes; j++)
            K(FS_NDOF*nodes[i]->index + 0, j) =  0.0;
        K(FS_NDOF*nodes[i]->index + 0, FS_NDOF*nodes[i]->index + 0) =  1.0;
    }
    if(nodes[i]->lockStatus[1])
    {
        for(int j=0; j<FS_NDOF*nNodes; j++)
            K(FS_NDOF*nodes[i]->index + 1, j) =  0.0;
        K(FS_NDOF*nodes[i]->index + 1, FS_NDOF*nodes[i]->index + 1) =  1.0;
    }
    if(nodes[i]->lockStatus[2])
    {
        for(int j=0; j<FS_NDOF*nNodes; j++)
            K(FS_NDOF*nodes[i]->index + 2, j) =  0.0;
        K(FS_NDOF*nodes[i]->index + 2, FS_NDOF*nodes[i]->index + 2) =  1.0;
    }
    if(nodes[i]->lockStatus[3])
    {
        for(int j=0; j<FS_NDOF*nNodes; j++)
            K(FS_NDOF*nodes[i]->index + 3, j) =  0.0;
        K(FS_NDOF*nodes[i]->index + 3, FS_NDOF*nodes[i]->index + 3) =  1.0;
    }
    if(nodes[i]->lockStatus[4])
    {
        for(int j=0; j<FS_NDOF*nNodes; j++)
            K(FS_NDOF*nodes[i]->index + 4, j) =  0.0;
        K(FS_NDOF*nodes[i]->index + 4, FS_NDOF*nodes[i]->index + 4) =  1.0;
    }
}


#pragma omp parallel for
for(int i=0; i<nNodes; i++)
{
    f(FS_NDOF*nodes[i]->index + 0) = nodes[i]->loadValues[0];
    f(FS_NDOF*nodes[i]->index + 1) = nodes[i]->loadValues[1];
    f(FS_NDOF*nodes[i]->index + 2) = nodes[i]->loadValues[2];
    f(FS_NDOF*nodes[i]->index + 3) = nodes[i]->loadValues[3];
    f(FS_NDOF*nodes[i]->index + 4) = nodes[i]->loadValues[4];
}


    std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;

    std::cout<<"\n\n dim = "<<K.n<<std::flush;

    K.solve(f, x);

    std::cout<<"\n\n *********** SOLVER LINEAR SYSTEM ***********"<<std::flush;
    w = new double[nNodes];
    for(int i=0; i<nNodes; i++)
        w[i] = x(5*i+0, 0);

    //std::cout<<"\n\n TEMPO = "<<double(clock() - t_start)/CLOCKS_PER_SEC;

    for(int i=0; i<nNodes; i++)
        std::cout<<"\n"<<i<<"\t"<<x(5*i+0, 0)<<"\t"<<x(5*i+1, 0)<<"\t"<<x(5*i+2, 0)<<"\t"<<x(5*i+3, 0)<<"\t"<<x(5*i+4, 0);
}
