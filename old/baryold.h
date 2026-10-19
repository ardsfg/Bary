#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <tuple>
using namespace std;

class Bary{
    private:
        double x,y,z;
        tuple<double,double,double> sides;
    public:
        Bary(double x, double y, double z, tuple<double,double,double> sides){
            this->x = x;
            this->y = y;
            this->z = z;
            this->sides = sides;
        }

        double getX() const{
            return x;
        }
        double getY()const {
            return y;
        }
        double getZ()const {
            return z;
        }
        tuple<double,double,double> getSides() const{
            return sides;
        }



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

        

};