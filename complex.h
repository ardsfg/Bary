#include <iostream>
#include <stdexcept>    
#include <vector>

#include <cmath>
#include <cstdlib>
#include <ctime>    
using namespace std;

class Complex{
private:
    double real;
    double imag;
public:
    Complex(double r, double i) : real(r), imag(i) {}
    double getReal() const { return real; }
    double getImag() const { return imag; }
    double getMag() const { return sqrt(real * real + imag * imag); }
    double getArg() const { return atan2(imag, real); }

Complex getConjugate(Complex a) { return Complex(a.getReal(), -a.getImag()); }

    Complex add(const Complex& a, const Complex& b) {
        return Complex(a.getReal() + b.getReal(), a.getImag() + b.getImag());
    }
    Complex subtract(const Complex& a, const Complex& b) {
        return Complex(a.getReal() - b.getReal(), a.getImag() - b.getImag());
    }
    Complex multiply(const Complex& a, const Complex& b) {
        return Complex(a.getReal() * b.getReal() - a.getImag() * b.getImag(),
                       a.getReal() * b.getImag() + a.getImag() * b.getReal());
    }
    Complex divide(const Complex& a, const Complex& b) {
        double sqmagnitude = b.getReal() * b.getReal() + b.getImag() * b.getImag();
        return Complex((a.getReal() * b.getReal() + a.getImag() * b.getImag())/sqmagnitude,
                       (a.getImag() * b.getReal() - a.getReal() * b.getImag())/sqmagnitude);
    }
Complex generateRandom(double r) {
    srand(time(0));
    double theta = (rand() / (double)RAND_MAX) * 2 * M_PI; // Random angle
    return Complex(r * cos(theta), r * sin(theta));
}

Complex getIntersection(const Complex& a, const Complex& b, const Complex& c, const Complex& d) {
        Complex x = subtract(a,b);
        Complex y = subtract(c,d);



        Complex wx=subtract(multiply(getConjugate(a),b),multiply(getConjugate(b),a));
        Complex wy=subtract(multiply(getConjugate(c),d),multiply(getConjugate(d),c));


        

        Complex num = subtract(multiply(wx,y), multiply(x,wy));
        Complex denom = subtract(multiply(getConjugate(x),y), multiply(getConjugate(y),x));

        Complex z = divide(num, denom);


        return z;
    }

Complex getMidpoint(const Complex& a, const Complex& b) {
        return Complex((a.getReal() + b.getReal()) / 2, (a.getImag() + b.getImag()) / 2);
    }


    Complex getReflection(Complex& a, Complex& b, const Complex& c) {
        Complex x=subtract(b,a);
        Complex y = subtract(multiply(a,getConjugate(b)), multiply(getConjugate(a),b));

        Complex z = divide(add(y,multiply(getConjugate(c),x)),getConjugate(x));

        return Complex(-z.getReal(), -z.getImag()); // Return the reflection point
    }



    bool isOrthogonal(Complex& a, Complex& b, Complex& c, Complex& d){
        Complex x=subtract(a,b);
        Complex y=subtract(c,d);

        if(divide(x, y).getArg() + getConjugate(divide(x, y)).getArg() < 1e-10) {
            return true; // The points are orthogonal
        }
        return false;
    }


    bool isCollinear(Complex& a, Complex& b, Complex& c) {
        Complex x = subtract(c, a);
        Complex y = subtract(c, b);


        Complex x1=getConjugate(x);
        Complex y1=getConjugate(y);

        Complex z = divide(divide(x,x1), divide(y,y1));

        if(abs(z.getArg())<1e-9) {
            return true; // The points are collinear
        }
        return false; // The points are not collinear
    }

    Complex rotate(Complex& a, Complex& b, double angle) {
        Complex x = subtract(b, a);

        return Complex(a.getReal()+x.getReal()*cos(angle)-x.getImag()*sin(angle), a.getImag()+x.getReal()*sin(angle)+x.getImag()*cos(angle));
    }

    Complex circumcenter(Complex& a, Complex& b, Complex& c) {


        Complex num1=add(multiply(multiply(add(a,b),getConjugate(subtract(b,a))),subtract(b,c)), multiply(subtract(b,c),subtract(multiply(getConjugate(a),b),multiply(getConjugate(b),a))));
        Complex num2=add(multiply(multiply(add(c,b),getConjugate(subtract(b,c))),subtract(b,a)), multiply(subtract(b,a),subtract(multiply(getConjugate(c),b),multiply(getConjugate(b),c))));
        Complex num = subtract(num1, num2);


        Complex denom1=multiply(getConjugate(subtract(b,a)),subtract(b,c));
        Complex denom2=multiply(getConjugate(subtract(b,c)),subtract(b,a));
        Complex denom = subtract(denom1, denom2);

        Complex d = divide(num, denom);

        return d;





    }


    vector<Complex> circleIntersection(Complex& o1, Complex& o2, double r1, double r2){
        vector<Complex> intersections;


        if(r1+r2 < abs(subtract(o1, o2).getMag())) {
            throw invalid_argument("The circles do not intersect.");
        }


        double d= subtract(o2, o1).getMag();
        double angle = acos((d*d+r1*r1-r2*r2)/(2*d*r1));

        Complex mid = subtract(o2, o1);
        mid = Complex(mid.getReal() / d, mid.getImag() / d);
        mid = Complex(mid.getReal() * r1, mid.getImag() * r1);
        mid = add(o1, mid);
    

        intersections.push_back(rotate(o1, mid, angle));
        intersections.push_back(rotate(o1, mid, -angle));

        return intersections;
    }

    bool isSame(Complex& a, Complex& b) {
        if(abs(a.getReal() - b.getReal()) < 1e-9 && abs(a.getImag() - b.getImag()) < 1e-9) {
            return true; // The points are the same
        }
        return false; // The points are not the same
    }

    bool isCyclic(Complex& a, Complex& b, Complex& c, Complex& d) {
        Complex x = subtract(c, a);
        Complex y = subtract(c, b);


        Complex x1 = subtract(d, a);
        Complex y1 = subtract(d, b);


        Complex z = divide(divide(x,x1), divide(y,y1));


        if(abs(z.getImag()/z.getReal()) < 1e9) {
            return true; // The points are cyclic
        }
        return false; // The points are not cyclic
    }


    Complex getFoot(Complex& a, Complex& b, Complex& c) {
        return getMidpoint(getReflection(a, b, c), c);
    }

    double getArea(Complex& a, Complex& b, Complex& c) {
        Complex x = subtract(b, a);
        Complex y = subtract(c, a);


        return 0.5 * abs(x.getReal() * y.getImag() - x.getImag() * y.getReal());
    }


};

