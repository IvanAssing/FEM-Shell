#ifndef POLYNOMIAL2D_H
#define POLYNOMIAL2D_H

#include "polynomial1d.h"
#include "functor2d.h"


    class Polynomial2D : public Functor2D
    {
        private:
            int n1;
            int n2;

            double *an;

            void resize(int degree1, int degree2);

        public:
            // constructors
            Polynomial2D();
            Polynomial2D(const Polynomial2D& obj);
            Polynomial2D(int degree1, int degree2, double *an);
            Polynomial2D(int degree, double *an);
            Polynomial2D(int degree1, int degree2, double alpha = 1.0);



            // operators
            virtual double operator()(double v1, double v2);

            Polynomial2D& operator*(const Polynomial2D& p1) const;
            Polynomial2D& operator=(const Polynomial2D& p1) ;
            Polynomial2D& operator+(const Polynomial2D& p1) const;
            Polynomial2D& operator-(const Polynomial2D& p1) const;
            //Rational2D& operator/(const Polynomial2D& obj) const;
            Polynomial2D& operator*(const double alpha) const;

            // methods

            int getDegree1(void) const;
            int getDegree2(void) const;

            double integral(double xinf = -1.0, double xsup = 1.0, double yinf = -1.0, double ysup = 1.0) const;
            Polynomial1D& integral1(double yinf = -1.0, double ysup = 1.0) const;
            Polynomial1D& integral2(double xinf = -1.0, double xsup = 1.0) const;

            Polynomial1D& evaluate1(double v1) const;
            Polynomial1D& evaluate2(double v2) const;

            Polynomial2D& differential1(int nDiffOrder = 1) const;
            Polynomial2D& differential2(int nDiffOrder = 1) const;

            friend std::ostream& operator<<(std::ostream& out, Polynomial2D& obj);
            friend Polynomial2D& operator*(double alpha, Polynomial2D& obj);

            // destructor
            virtual ~Polynomial2D();
    };

#endif // POLYNOMIAL2D_H
