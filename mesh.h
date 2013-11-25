#ifndef MESH_H
#define MESH_H

#include "node.h"
#include "elementdkt.h"


class Mesh
{
    public:

        Node **nodes;
        ElementDKT **elements;

        double *w;

        int nNodes;
        int nElements;

        void plot(Matrix &f);

        void draw(void);

        Mesh();
};

#endif // MESH_H
