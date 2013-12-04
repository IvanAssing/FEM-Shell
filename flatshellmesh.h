#ifndef FLATSHELLMESH_H
#define FLATSHELLMESH_H

#include "mesh.h"
#include "node.h"
#include "elementsqn.h"
#include "lagrange.h"
#include "polynomial2d.h"

class FlatShellMesh : public Mesh
{
    public:
        Lagrange *L;
        Polynomial2D **BftDBf;
        Polynomial2D **BctBc;
        Polynomial2D **BmtDBm;

        int npx, npy;


        FlatShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSQN **_elements, int _npx, int _npy, Matrix _Df, double _GKt, Matrix _Dm);

        Node **nodes;
        ElementSQN **elements;
        Matrix Df, Dm;
        double GKt;

        double *w;

        int nNodes;
        int nElements;

        virtual void draw(void);

    private:
        int npt;
};

#endif // FLATSHELLMESH_H
