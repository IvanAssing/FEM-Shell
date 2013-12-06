#include "polynomial2d.h"


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>


Polynomial2D::Polynomial2D()
    :n1(1),n2(1)
{
    an = new double[1];
    an[0] = 0.0;
}


Polynomial2D::Polynomial2D(const Polynomial2D& obj)
    :n1(obj.n1), n2(obj.n2)
{
    an = new double[n1*n2];
    for(int i=0; i<n1*n2; i++)
        an[i] = obj.an[i];
}


Polynomial2D::Polynomial2D(int degree1, int degree2, double*coefficients)
    :n1(degree1+1), n2(degree2+1)
{
    an = new double[n1*n2];

//#pragma omp parallel for
    for(int i=0; i<n1*n2; i++)
        an[i] = coefficients[i];
}


Polynomial2D::Polynomial2D(int degree1, int degree2, double alpha)
    :n1(degree1+1), n2(degree2+1)
{
    an = new double[n1*n2];
    for(int i=0; i<n1*n2; i++)
        an[i] = alpha;
}


Polynomial2D::Polynomial2D(int degree, double*coefficientes)
    :n1(degree+1), n2(degree+1)
{
    an = new double[n1*n2];

//#pragma omp parallel for
    for(int i=0; i<n1*n2; i++)
        an[i] = coefficientes[i];
}


Polynomial2D::~Polynomial2D()
{
    //delete [] an;
}


double Polynomial2D::operator()(double v1, double v2)
{
    double*ptr = &an[0];
    double sum = *(ptr++);

    double xn,yn;

    xn = 1.0;
    for(int i=1; i<n1; i++)
    {
        xn *= v1;
        sum += *(ptr++)*xn;
    }

    yn = 1.0;
    for(int i=1; i<n2; i++)
    {
        yn *= v2;
        xn = 1.0;
        sum += *(ptr++)*yn;
        for(int j=1; j<n1; j++)
        {
            xn *= v1;
            sum += *(ptr++)*xn*yn;
        }
    }

    return sum;
}


Polynomial2D& Polynomial2D::operator+(const Polynomial2D& obj) const
{
    Polynomial2D *sum = new Polynomial2D;

    sum->resize(std::max(n1,obj.n1), std::max(n2,obj.n2));

    // THIS
//#pragma omp parallel for
    for(int i=0; i<n2; i++)
        for(int j=0; j<n1; j++)
            sum->an[i*sum->n1+j] += an[i*n1+j];

    // OBJ
//#pragma omp parallel for
    for(int i=0; i<obj.n2; i++)
        for(int j=0; j<obj.n1; j++)
            sum->an[i*sum->n1+j] += obj.an[i*obj.n1+j];

    return *sum;

}


Polynomial2D& Polynomial2D::operator-(const Polynomial2D& obj) const
{
    Polynomial2D *sum = new Polynomial2D;

    sum->resize(std::max(n1,obj.n1), std::max(n2,obj.n2));

    // THIS
    for(int i=0; i<n2; i++)
        for(int j=0; j<n1; j++)
            sum->an[i*sum->n1+j] += an[i*n1+j];
    // OBJ
    for(int i=0; i<obj.n2; i++)
        for(int j=0; j<obj.n1; j++)
            sum->an[i*sum->n1+j] -= obj.an[i*obj.n1+j];

    return *sum;

}


Polynomial2D& Polynomial2D::operator*(const double alpha) const
{
    Polynomial2D *prod = new Polynomial2D(*this);

    for(int i=0; i<prod->n1*n2; i++)
        prod->an[i] *= alpha;

    return *prod;
}


Polynomial2D& Polynomial2D::operator*(const Polynomial2D& obj) const
{
    Polynomial2D *prod = new Polynomial2D;

    prod->resize(n1 + obj.n1 - 1, n2 + obj.n2 - 1);

    for(int it=0; it<n2; it++)
        for(int jt=0; jt<n1; jt++)
            for(int i=0; i<obj.n2; i++)
                for(int j=0; j<obj.n1; j++)
                    prod->an[it*prod->n1+jt + i*prod->n1+j] += an[it*n1+jt] * obj.an[i*obj.n1+j];

    return *prod;
}


Polynomial2D& Polynomial2D::operator=(const Polynomial2D& obj)
{
    this->resize(obj.n1, obj.n2);

    for(int i=0; i<obj.n1*obj.n2; i++)
        an[i] = obj.an[i];

    return *this;
}


//Rational2D& Polynomial2D::operator/(const Polynomial2D& obj) const
//{
//    Rational2D *ratio = new Rational2D(*this, obj);

//    return *ratio;
//}


void Polynomial2D::resize(int newSize1, int newSize2)
{
    this->n1 = newSize1;
    this->n2 = newSize2;

    delete [] an;

    an = new double[n1*n2];

    for(int i=0; i<n1*n2; i++)
        an[i] = 0.0;

}


int Polynomial2D::getDegree1(void) const
{
    return n1-1;
}


int Polynomial2D::getDegree2(void) const
{
    return n2-1;
}


Polynomial1D& Polynomial2D::integral1(double xinf, double xsup) const
{
    double*intAn = new double[n2];

    double dn;
    double sum;
    double xa;
    double xb;


    for(int k=0; k<n2; k++){
        xa = xinf;
        xb = xsup;
        dn = 1.0;
        sum = an[k*n2]*(xb-xa);

        for(int i=1; i<n1; i++){
            xa *= xinf;
            xb *= xsup;
            dn += 1.0;
            sum += an[k*n2+i]*(xb-xa)/dn;
        }

        intAn[k] = sum;
    }

    Polynomial1D *intX = new Polynomial1D(n2-1, intAn);

    delete [] intAn;

    return *intX;
}


Polynomial1D& Polynomial2D::integral2(double yinf, double ysup) const
{
    double*intAn = new double[n1];

    double dn;
    double sum;
    double ya;
    double yb;


    for(int k=0; k<n1; k++){
        ya = yinf;
        yb = ysup;
        dn = 1.0;
        sum = an[k]*(yb-ya);

        for(int i=1; i<n2; i++){
            ya *= yinf;
            yb *= ysup;
            dn += 1.0;
            sum += an[k+i*n2]*(yb-ya)/dn;
        }

        intAn[k] = sum;
    }
    Polynomial1D *intY = new Polynomial1D(n1-1, intAn);

    delete [] intAn;

    return *intY;
}


double Polynomial2D::integral(double xinf, double xsup, double yinf, double ysup) const
{
    double*intAn = new double[n2];

    double dn;
    double sum;
    double xa;
    double xb;


    for(int k=0; k<n2; k++){
        xa = xinf;
        xb = xsup;
        dn = 1.0;
        sum = an[k*n2]*(xb-xa);

        for(int i=1; i<n1; i++){
            xa *= xinf;
            xb *= xsup;
            dn += 1.0;
            sum += an[k*n2+i]*(xb-xa)/dn;
        }

        intAn[k] = sum;
    }
    Polynomial1D *intX = new Polynomial1D(n2-1, intAn);

    delete [] intAn;

    return intX->integral(yinf,ysup);
}


Polynomial2D& Polynomial2D::differential1(int nDiffOrder) const
{
    if(nDiffOrder<1) return *new Polynomial2D(*this);

    Polynomial2D *out = new Polynomial2D;

    if(nDiffOrder>=n1) return *out;

    out->resize(n1 - nDiffOrder, n2);

    for(int k=0; k<n2; k++){
        double fat = 1.0;

        for(int i=0; i<nDiffOrder; i++)
            fat *= static_cast<double>(n1 - i - 1);

        for(int i=n1-1; i-nDiffOrder>=0; i--){
            out->an[k*out->n1+i-nDiffOrder] = an[k*n2+i]*fat;
            fat *= static_cast<double>(i-nDiffOrder)/static_cast<double>(i);
        }
    }

    return *out;
}


Polynomial2D& Polynomial2D::differential2(int nDiffOrder) const
{
    if(nDiffOrder<1) return *new Polynomial2D(*this);

    Polynomial2D *out = new Polynomial2D;

    if(nDiffOrder>=n2) return *out;

    out->resize(n1, n2 - nDiffOrder);

    for(int k=0; k<n1; k++){
        double fat = 1.0;

        for(int i=0; i<nDiffOrder; i++)
            fat *= static_cast<double>(n2 - i - 1);

        for(int i=n2-1; i-nDiffOrder>=0; i--){
            out->an[k+(i-nDiffOrder)*n2] = an[k+i*n2]*fat;
            fat *= static_cast<double>(i-nDiffOrder)/static_cast<double>(i);
        }
    }

    return *out;
}


Polynomial1D& Polynomial2D::evaluate1(double x) const
{
    double *an = new double[n2];


    for(int i=0; i<n2; i++){
        an[i] = this->an[i*n1];
        double fatx = x;
        for(int j=1; j<n1; j++){
            an[i] += an[i*n1 + j]*fatx;
            fatx *= x;
        }
    }

    Polynomial1D *obj = new Polynomial1D(n2-1, an);

    delete [] an;

    return *obj;
}


Polynomial1D& Polynomial2D::evaluate2(double y) const
{
    double *an = new double[n1];


    for(int j=0; j<n1; j++)
        an[j] = an[j];

    double faty = y;
    for(int i=1; i<n2; i++){
        for(int j=0; j<n1; j++)
            an[j] += an[j + i*n1]*faty;
        faty *= y;
    }

    Polynomial1D *obj = new Polynomial1D(n1-1, an);

    delete [] an;

    return *obj;
}



std::ostream&  operator<<(std::ostream& out, Polynomial2D& obj)
{
    out<<std::setprecision(4)<<std::setw(10)<<std::uppercase;
    out<<std::scientific<<obj.an[0];
    for(int i=1; i<obj.n1; i++)
        if(fabs(obj.an[i])>1.0e-5)
            out<<" + "<<obj.an[i]<<" X^"<<i;
    for(int i=1; i<obj.n2; i++){
        if(fabs(obj.an[i*obj.n1])>1.0e-5)
            out<<" + "<<obj.an[i*obj.n1]<<" Y^"<<i;
        for(int j=1; j<obj.n1; j++){
            if(fabs(obj.an[i*obj.n1+j])>1.0e-5)
                out<<" + "<<obj.an[i*obj.n1+j]<<" X^"<<j<<" Y^"<<i;
        }
    }

    return out;
}

Polynomial2D& operator*(double alpha, Polynomial2D& obj)
{
    return obj*alpha;
}

