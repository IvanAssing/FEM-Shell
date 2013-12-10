#ifndef THINSHELLMESH_H
#define THINSHELLMESH_H


#include "mesh.h"
#include "node.h"
#include "elementdkt.h"
#include "matrix.h"


class ThinShellMesh : public Mesh
{
    public:
        ThinShellMesh(){}
        ThinShellMesh(int _nNodes, Node ** _nodes, int _nElements, ElementDKT **_elements, Matrix _D);

    public:

        Node **nodes;
        ElementDKT **elements;
        Matrix D;

        int nNodes;
        int nElements;

        double **results;


        void plot(void);

        virtual void draw(DataGraphic &data);
        virtual void solve(void);

};
#endif // THINSHELLMESH_H
