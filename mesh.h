#ifndef MESH_H
#define MESH_H

#include "node.h"
#include "element.h"


enum vout{
    U,
    V,
    W,
    RX,
    RY,
    RZ,
    FX,
    FY,
    FZ,
    MX,
    MY,
    MZ
};

class Mesh
{
    public:
        int nNodes;
        int nElements;

        Node **nodes;
        Element **elements;

    public:
        Mesh(){}
        virtual void plot(void){}
        virtual void draw(vout option)=0;
        virtual void solve(void)=0;
};


void getMaxMin(double *vector, int size, double &max, double &min);

#endif // MESH_H
