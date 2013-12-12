#ifndef SHELLMESH_H
#define SHELLMESH_H

#include "mesh.h"
#include "node.h"
#include "elementssqn.h"
#include "lagrange.h"
#include "polynomial2d.h"

class ShellMesh : public Mesh
{
    public:
        Lagrange *L;
        Polynomial2D **BftDBf;
        Polynomial2D **BctBc;
        Polynomial2D **BmtDBm;
        int npx, npy;

        Node **nodes;
        ElementSSQN **elements;
        Matrix Df, Dm;
        double GKt;
        double **results;
        int nNodes;
        int nElements;


        ShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementSSQN **_elements, int _npx, int _npy, Matrix _Df, double _GKt, Matrix _Dm);

        virtual void draw(DataGraphic &data);
        virtual void solve(void);

    private:
        int npt;
};

#endif // SHELLMESH_H
