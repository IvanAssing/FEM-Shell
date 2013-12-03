#ifndef MESH_H
#define MESH_H

#include "node.h"
#include "element.h"


class Mesh
{
    public:
        int nNodes;
        int nElements;

        Node **nodes;
        Element **elements;

    public:
        Mesh();
        virtual void plot(void){}
        virtual void draw(void)=0;
};

#endif // MESH_H
