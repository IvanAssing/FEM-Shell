#ifndef RATIONAL2D_H
#define RATIONAL2D_H

#include <iostream>

#include "functor2d.h"
#include "polynomial2d.h"


class Rational2D : public Functor2D
{
    private:
        Polynomial2D pnum;
        Polynomial2D pden;



    public:
        // constructors
        Rational2D();
        Rational2D(const Rational2D& obj);
        Rational2D(const Polynomial2D& numerator, const Polynomial2D& denominator);


        // operators
        virtual double operator()(double v1, double v2);

        Rational2D& operator*(const Rational2D& obj) const;
        Rational2D& operator=(const Rational2D& obj) ;
        Rational2D& operator+(const Rational2D& obj) const;
        Rational2D& operator-(const Rational2D& obj) const;
        Rational2D& operator*(const double alpha) const;

        // methods

        //INLINE int getOrderX(void) const;
        //INLINE int getOrderY(void) const;

        //real evalIntegral(real xinf = -1.0, real xsup = 1.0, real yinf = -1.0, real ysup = 1.0) const;
        //FXPolynomial& evalIntegralY(real yinf = -1.0, real ysup = 1.0) const;
        //FXPolynomial& evalIntegralX(real xinf = -1.0, real xsup = 1.0) const;

        //Polynomial2D& diffX(int nDiffOrder = 1) const;
        //Polynomial2D& diffY(int nDiffOrder = 1) const;

        friend std::ostream& operator<<(std::ostream& out, Rational2D& obj);

        // destructor
        virtual ~Rational2D();
};


#endif // RATIONAL2D_H
