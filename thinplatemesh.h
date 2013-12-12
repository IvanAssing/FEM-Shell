#ifndef THINPLATEMESH_H
#define THINPLATEMESH_H

#include "mesh.h"
#include "node.h"
#include "elementdkt.h"
#include "matrix.h"


class ThinPlateMesh : public Mesh
{
    public:
        Node **nodes;
        ElementDKT **elements;
        Matrix D;

        int nNodes;
        int nElements;
        double **results;
        Gnuplot *gnuplot;

        ThinPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementDKT **_elements, Matrix _D);

        void plot(vout data);
        virtual void draw(DataGraphic &data);
        virtual void solve(void);
};

#endif // THINPLATEMESH_H
