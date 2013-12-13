#ifndef POLYNOMIAL1D_H
#define POLYNOMIAL1D_H

#include <iostream>



class Polynomial1D
{
        //private:
    public:
        double *an;
        int n;

        void resize(int n);

    public:
        // constructors
        Polynomial1D();
        Polynomial1D(const Polynomial1D& b);
        Polynomial1D(int order, double *an);

        // operators
        double operator()(double x) const;
        Polynomial1D& operator*(const Polynomial1D& b) const;
        Polynomial1D& operator=(const Polynomial1D& b) ;
        Polynomial1D& operator+(const Polynomial1D& b) const;
        Polynomial1D& operator-(const Polynomial1D& b) const;
        Polynomial1D& operator*(const double alpha) const;

        // methods
        int getDegree(void) const;
        double integral(double xinf = -1.0, double xsup = 1.0) const;
        Polynomial1D& differential(int nDiffOrder = 1) const;

        friend std::ostream& operator<<(std::ostream& out, Polynomial1D& b);

        // destructor
        virtual ~Polynomial1D();
};


#endif // POLYNOMIAL1D_H
