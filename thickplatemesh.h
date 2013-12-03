#ifndef THICKPLATEMESH_H
#define THICKPLATEMESH_H

#include "node.h"
#include "elementqn.h"
#include "lagrange.h"
#include "polynomial2d.h"

class ThickPlateMesh
{
    public:
        Lagrange *L;
        Polynomial2D **BftDBf;
        Polynomial2D **BctBc;

        int npx, npy;


        ThickPlateMesh(int _nNodes, Node ** _nodes, int _nElements, ElementQN **_elements, int _npx, int _npy);

        Node **nodes;
        ElementQN **elements;

        double *w;

        int nNodes;
        int nElements;

        void draw(void);

    private:
        int npt;
};

#endif // THICKPLATEMESH_H
