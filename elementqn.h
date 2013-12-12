#ifndef ELEMENTQN_H
#define ELEMENTQN_H

#include "lagrange.h"
#include "node.h"
#include "matrix.h"
#include "element.h"

class ElementQN : public Element
{
    public:
        int np;
        int npx, npy;
        bool selectiveIntegration;
        Polynomial2D J;
        Node **nodes;

        ElementQN(int npx_, int npy_, Node **nodes_, bool _selectiveIntegracion = false);

        void getStiffnessMatrix(Matrix &k, Polynomial2D **Bf, Polynomial2D **Bc, Lagrange *L);
        void evalResults(Matrix &M, Matrix &Q, Matrix &U, Polynomial2D **Bf, Polynomial2D **Bc);

        virtual void draw(void);
};

#endif // ELEMENTQN_H
