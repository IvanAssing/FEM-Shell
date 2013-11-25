#include "mesh2.h"

Mesh2::Mesh2()
{

    int np = 3+1;
    int ne = 5;
    int nx = (ne-1)*np, ny = (ne-1)*np;

    nNodes = nx*ny;

    nodes = new Node*[nNodes];

    double dx = 1.0/(nx-1);
    double dy = 1.0/(ny-1);


    for(int i=0; i<ny; i++)
        for(int j=0; j<nx; j++)
            nodes[j+nx*i] = new Node(j+nx*i, j*dx, i*dy);

    Node **ptrNodes = new Node*[np*np];


    nElements = ne*ne;
    elements = new ElementQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ne; ie++)
        for(int je=0; je<ne; je++)
        {
            ni = 0;
            for(int i=0; i<np; i++)
                for(int j=0; j<np; j++)
                    ptrNodes[ni++] = nodes[ie*(np-1)*nx + i*nx + je*(np-1) + j];
            elements[elementIndex++] = new ElementQN(ptrNodes);
        }

    Lagrange L(3, 3);

    std::cout<<L.N[15](1, 1);

    Polynomial2D Bf[3][3*16];
    Polynomial2D Bc[2][3*16];

    for(int i=0; i<16; i++)
    {
        Bf[0][3*i+2] = -1.0*L.D1[i];
        Bf[1][3*i+1] = L.D2[i];
        Bf[2][3*i+1] = L.D1[i];
        Bf[2][3*i+2] = -1.0*L.D2[i];

        Bc[0][3*i] = L.D1[i];
        Bc[0][3*i+2] = L.N[i];
        Bc[1][3*i] = L.D2[i];
        Bc[1][3*i+1] = -1.0*L.N[i];
    }



    double vi = 0.3;
    double E = 200e9;
    double t = 0.02;

    double GKt = 10.;

    double Ept = E*t*t*t/(12.*(1.0-vi*vi));

    Matrix D(3,3);

    D(0, 0, Ept);
    D(0, 1, Ept*vi);
    D(1, 0, Ept*vi);
    D(1, 1, Ept);
    D(2, 2, Ept*(1-vi)/2.0);



    Polynomial2D BftDBf[3*16][3*16];
    Polynomial2D DB[3][3*16];

    for(int i=0; i<3; i++)
        for(int j=0; j<3*16; j++)
            DB[i][j] = Bf[0][j]*D(i,0) + Bf[1][j]*D(i,1) + Bf[2][j]*D(i,2);

    for(int i=0; i<3*16; i++)
        for(int j=0; j<3*16; j++)
            BftDBf[i][j] = Bf[0][i]*DB[0][j] + Bf[1][i]*DB[1][j] + Bf[2][i]*DB[2][j];

    Polynomial2D BctBc[3*16][3*16];

    for(int i=0; i<3*16; i++)
        for(int j=0; j<3*16; j++)
            BctBc[i][j] = (Bc[0][i]*Bc[0][j] + Bc[1][i]*Bc[1][j])*GKt;


    Matrix K(3*nNodes, 3*nNodes);

    Matrix f(3*nNodes, 1);
    Matrix x(3*nNodes, 1);

    elements[0]->getStiffnessMatrix(K, &BftDBf, &BctBc);






}
