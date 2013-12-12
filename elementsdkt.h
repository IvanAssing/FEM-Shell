#ifndef ELEMENTSDKT_H
#define ELEMENTSDKT_H

#include "node.h"
#include "polynomial2d.h"
#include <iostream>
#include "matrix.h"
#include "element.h"

class ElementSDKT : public Element
{
    public:
        int index;
        Node *n1, *n2, *n3;
        Polynomial2D **Bf;
        double **Bm;

        ElementSDKT(int index, Node *node1, Node *node2, Node *node3);

        void evaluateTransformationMatrix(void);
        void getStiffnessMatrix(Matrix &k, Matrix &Df, Matrix &Dm);

        void evalResults(Matrix &M, Matrix &U, Matrix &D);

        virtual void draw(void);
};


#endif // ELEMENTSDKT_H
