#include <iostream>
#include "complex.h"
#include "bary.h"
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <iomanip>
#include <cmath>
using namespace std;


double determinant(Bary A, Bary B, Bary C){
            double ret = abs(A.getX()*(B.getY()*C.getZ() - C.getY()*B.getZ()) - A.getY()*(B.getX()*C.getZ() - C.getX()*B.getZ()) + A.getZ()*(B.getX()*C.getY() - C.getX()*B.getY()));
            ret/=(A.getX()+A.getY()+A.getZ())*(B.getX()+B.getY()+B.getZ())*(C.getX()+C.getY()+C.getZ());
            return ret;
        }
        bool isCollinearBary(Bary A, Bary B, Bary C){
            if(determinant(A, B, C) < 1e-9){
                return true;
            }
            return false;
        }

        Bary isogonalConjugate(Bary A){
            Bary ret = Bary(get<0>(A.getSides())*get<0>(A.getSides())/A.getX(), get<1>(A.getSides())*get<1>(A.getSides())/A.getY(), get<2>(A.getSides())*get<2>(A.getSides())/A.getZ(), A.getSides());
            double x = ret.getX()/(ret.getX() + ret.getY() + ret.getZ());
            double y = ret.getY()/(ret.getX() + ret.getY() + ret.getZ());
            double z = ret.getZ()/(ret.getX() + ret.getY() + ret.getZ());
            ret = Bary(x, y, z, A.getSides());
            return ret;
        }
        Bary isotomicConjugate(Bary A){
            Bary ret = Bary(1.0/get<0>(A.getSides()), 1.0/get<1>(A.getSides()), 1.0/get<2>(A.getSides()), A.getSides());
            double x = ret.getX()/(ret.getX() + ret.getY() + ret.getZ());
            double y = ret.getY()/(ret.getX() + ret.getY() + ret.getZ());
            double z = ret.getZ()/(ret.getX() + ret.getY() + ret.getZ());
            ret = Bary(x, y, z, A.getSides());
            return ret;
        }

        //circumcenter = (a^2(b^2+c^2-a^2) : b^2(c^2+a^2-b^2) : c^2(a^2+b^2-c^2))
        //orthocenter = (S_b*S_c : S_c*S_a : S_a*S_b)



        Bary intersection(Bary A, tuple<double,double,double> t1, tuple<double,double,double> t2){
            // Normalize the input barycentric coordinates (use 't1' tuple as the point); return as a Bary with the current triangle sides.
            if((get<1>(t1) - get<0>(t1))/(get<1>(t2)-get<0>(t2)) - (get<2>(t1) - get<0>(t1))/(get<2>(t2)-get<0>(t2)) < 1e-12){
                throw invalid_argument("The lines are parallel and do not intersect.");
            }


            double w1 = get<0>(t1);
            double w2 = get<1>(t1);
            double w3 = get<2>(t1);
            double v1 = get<0>(t2);
            double v2 = get<1>(t2);
            double v3 = get<2>(t2);


            double z = ((v1*(w2-w1)-w1)/(v2-v1))/(w3-w1-(v3-v1)*(w2-w1)/(v2-v1));
            double y = (-v1*(w2-w1)/(v2-v1) - z*(v3-v1)*(w2-w1)/(v2-v1))/(w2-w1);
            double x = 1 - y - z;

            return Bary(x, y, z, A.getSides());



        }

        double distance(Bary A, Bary B){
            double x = (get<0>(A.getSides())*(A.getX()/(A.getX() + A.getY() + A.getZ()) - B.getX()/(B.getX() + B.getY() + B.getZ())));
            double y = (get<1>(A.getSides())*(A.getY()/(A.getX() + A.getY() + A.getZ()) - B.getY()/(B.getX() + B.getY() + B.getZ())));
            double z = (get<2>(A.getSides())*(A.getZ()/(A.getX() + A.getY() + A.getZ()) - B.getZ()/(B.getX() + B.getY() + B.getZ())));
            return sqrt(x*x + y*y + z*z);
        }

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





int main(){
     // Seed the random number generator
     /* complex stuff apmo i think
    Complex o = Complex(0, 0);
    
    Complex a = Complex(0.98931871426575135331, 0.14576858922127425289);
    std::srand(std::time(0));
    Complex b = Complex(-0.76827002975782410221, 0.64012589494248095523);
    Complex c = Complex(0.17253415239780264812, -0.9850035361643996179);

    Complex h = add(a, add(b, c));
    Complex d = getReflection(a,c,h);
    Complex e = getReflection(a,b,h);

    std::cout << std::setprecision(20) << "a= (" << a.getReal() << ", " << a.getImag() << ")" << std::endl;
    std::cout << "b= (" << b.getReal() << ", " << b.getImag() << ")" << std::endl;
    std::cout << "c= (" << c.getReal() << ", " << c.getImag() << ")" << std::endl;

    d=Complex(-d.getReal(), -d.getImag());
    e=Complex(-e.getReal(), -e.getImag());

    std::cout<< "d= (" << d.getReal() << ", " << d.getImag() << ")" << std::endl;
    std::cout<< "e= (" << e.getReal() << ", " << e.getImag() << ")" << std::endl;




    Complex x = getIntersection(a,b,c,o);
    Complex y = getIntersection(a,c,b,o);

    std::cout << "x= (" << x.getReal() << ", " << x.getImag() << ")" << std::endl;
    std::cout << "y= (" << y.getReal() << ", " << y.getImag() << ")" << std::endl;   


    Complex o1=circumcenter(b,x,d);
    Complex o2=circumcenter(c,y,e);


    std::cout << "o1= (" << o1.getReal() << ", " << o1.getImag() << ")" << std::endl;
    std::cout << "o2= (" << o2.getReal() << ", " << o2.getImag() << ")" << std::endl; 


    std::vector<Complex> intersections = circleIntersection(o1, o2, subtract(o1,b).getMag(), subtract(o2,c).getMag());

    bool collinear=false;
    for(std::vector<Complex>::iterator it = intersections.begin(); it != intersections.end(); ++it) {

        Complex p= *it;

        std::cout << "Intersection point p= (" << p.getReal() << ", " << p.getImag() << ")" << std::endl;
        if(isCollinear(p, a, o)) {
            collinear = true;
        }
    }

    if(collinear) {
        std::cout << "The points are collinear." << std::endl;
    } else {
        std::cout << "The points are not collinear." << std::endl;
    }



    */
    tuple<double,double,double> side = make_tuple(3.0,4.0,5.0);
    cout<<determinant(Bary(1,2,3,side),Bary(3,1,2,side), Bary(2,3,1,side))<<endl;
    if(isCollinearBary(Bary(1,2,3,side),Bary(3,4,6,side), Bary(2,5,7.5,side))){
        cout<<"collinear"<<endl;
    } else {
        cout<<"not collinear"<<endl;
    }

}