#ifndef MESH_H
#define MESH_H

#include "node.h"
#include "element.h"

#include "gnuplot.h"


enum vout{
    U,
    V,
    W,
    RX,
    RY,
    RZ,
    MX,
    MY,
    MXY,
    QX,
    QY
};

struct DataGraphic{
        vout var;
        bool nodes;
        bool elements;
        bool load;
        bool def;
        double factor;
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
        virtual void plot(vout){}
        virtual void draw(DataGraphic &data)=0;
        virtual void solve(void)=0;
};


void getMaxMin(double *vector, int size, double &max, double &min);

#endif // MESH_H
