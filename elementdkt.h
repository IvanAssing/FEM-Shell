#ifndef ELEMENTDKT_H
#define ELEMENTDKT_H

#include "node.h"
#include "polynomial2d.h"
#include <iostream>
#include "matrix.h"

class ElementDKT
{
    public:
        int index;
        Node *n1, *n2, *n3;
        Polynomial2D **B;
        ElementDKT(int index, Node *node1, Node *node2, Node *node3);



        void evaluateTransformationMatrix(void);
        void getStiffnessMatrix(Matrix &k, Matrix &D);

        void draw(void);
};

#endif // ELEMENTDKT_H
