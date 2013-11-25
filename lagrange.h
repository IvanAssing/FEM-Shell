#ifndef LAGRANGE_H
#define LAGRANGE_H


#include "polynomial2d.h"


class Lagrange
{
    public:
        int n1, n2;
        int ne;
        Polynomial2D *N;
        Polynomial2D *D1;
        Polynomial2D *D2;



        Lagrange(int n1, int n2);


};

#endif // LAGRANGE_H
