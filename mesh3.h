#ifndef MESH3_H
#define MESH3_H

#include "matrix.h"
#include "node.h"
#include "elementsqn.h"

class Mesh3
{
    public:
        Mesh3();


        Node **nodes;
        ElementSQN **elements;

        double *w;

        int nNodes;
        int nElements;

        void plot(Matrix &f);

        void draw(void);

};

#endif // MESH2_H
