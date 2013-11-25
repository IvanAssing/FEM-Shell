#ifndef MESH2_H
#define MESH2_H

#include "matrix.h"
#include "node.h"
#include "elementqn.h"

class Mesh2
{
    public:
        Mesh2();


        Node **nodes;
        ElementQN **elements;

        double *w;

        int nNodes;
        int nElements;

};

#endif // MESH2_H