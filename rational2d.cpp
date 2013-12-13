#include "rational2d.h"

Rational2D::Rational2D()
{
    pnum = Polynomial2D(0,0,0.0);
    pden = Polynomial2D(0,0);
}


Rational2D::Rational2D(const Rational2D& obj)
{
    pnum = obj.pnum;
    pden = obj.pden;
}



Rational2D::Rational2D(const Polynomial2D& numerator, const Polynomial2D& denominator)
{
    pnum = numerator;
    pden = denominator;
}


double Rational2D::operator()(double v1, double v2)
{
    return pnum(v1,v2)/pden(v1,v2);
}


Rational2D& Rational2D::operator*(const Rational2D& obj) const
{
    Rational2D *prod = new Rational2D(pnum*obj.pnum, pden*obj.pden);
    return *prod;
}


Rational2D& Rational2D::operator+(const Rational2D& obj) const
{
    Rational2D *sum = new Rational2D;

    sum->pnum = pnum*obj.pden + pden*obj.pnum;
    sum->pden = pden*obj.pden;

    return *sum;
}


Rational2D& Rational2D::operator-(const Rational2D& obj) const
{
    Rational2D *sub = new Rational2D;

    sub->pnum = pnum*obj.pden - pden*obj.pnum;
    sub->pden = pden*obj.pden;

    return *sub;
}

Rational2D& Rational2D::operator=(const Rational2D& obj)
{
    this->pnum = obj.pnum;
    this->pden = obj.pden;

    return *this;
}

Rational2D& Rational2D::operator*(const double alpha) const
{
    Rational2D *prod = new Rational2D(*this);

    prod->pnum = pnum*alpha;

    return *prod;
}


Rational2D::~Rational2D()
{

}



std::ostream& operator<<(std::ostream& out, Rational2D& obj)
{
    out<<"( "<<obj.pnum<<" / "<<obj.pden<<" )";
    return out;
}




