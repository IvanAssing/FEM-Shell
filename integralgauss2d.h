#ifndef INTEGRALGAUSS2D_H
#define INTEGRALGAUSS2D_H

#include "functor2d.h"

extern "C" {
#include "quadmath.h" // Biblioteca para Precisão Quadrúpla
}

class IntegralGauss2D
{
    public:
        int n;

        IntegralGauss2D();
        IntegralGauss2D(int n);

        double operator()(Functor2D &f);

        static double int1P(Functor2D &function);
        static double int2P(Functor2D &function);
        static double int3P(Functor2D &function);
        static double int4P(Functor2D &function);
        static double int5P(Functor2D &function);
        static double int6P(Functor2D &function);
        static double int7P(Functor2D &function);
        static double int8P(Functor2D &function);
        static double int9P(Functor2D &function);
        static double int10P(Functor2D &function);

        static double intNP(int nPoints, Functor2D &function);

    private:
        __float128 *p;
        __float128 *w;

    public:
        virtual ~IntegralGauss2D();

};

#endif // INTEGRALGAUSS2D_H
