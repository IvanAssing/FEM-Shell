#ifndef THINPLATEMESH_H
#define THINPLATEMESH_H

#include "mesh.h"
#include "node.h"
#include "elementdkt.h"
#include "matrix.h"


class ThinPlateMesh : public Mesh
{
    public:
        ThinPlateMesh(){}
        ThinPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementDKT **_elements, Matrix _D);

    public:

        Node **nodes;
        ElementDKT **elements;
        Matrix D;

        int nNodes;
        int nElements;

        double **results;


        void plot(Matrix &f);

        virtual void draw(vout option);
        virtual void solve(void);

};

#endif // THINPLATEMESH_H
