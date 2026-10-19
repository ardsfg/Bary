#include "bits.stdc++.h"

class Bary{
    private:
        double a, b, c; 
        vector<vector<vector<int>>> numerator, denominator; // apow, bpow, cpow
    public: 
        Bary(double x, double y, double z){
            a=x;
            b=y;
            c=z;
            for(int i = 0; i<10; i++){
                vector<vector<int>> v;
                numerator.push_back(v);
                denominator.push_back(v);
                for(int j = 0; j<10; j++){
                    vector<int> v2;
                    numerator[i].push_back(v2);
                    denominator[i].push_back(v2);
                }
            }
        }
        double getA(){
            return a;
        }
        double getB(){
            return b;
        }
        double getC(){
            return c;
        }
        tuple<poly,poly,poly> normalize(tuple<poly,poly,poly> p){
            
        }
        tuple<poly,poly,poly> add(tuple<poly,poly,poly> p1, tuple<poly,poly,poly> p2){

        }
        tuple<poly,poly,poly> dilate(tuple<poly,poly,poly> p1, tuple<poly,poly,poly> p2, double d){
            return multiply(normalize(p2)-normalize(p1), d) + p1;
        }
        tuple<poly, poly, poly> perpendicular_line(tuple<poly,poly,poly> line, tuple<poly,poly,poly> point){
            
        }
        tuple<poly,poly,poly> line(tuple<poly,poly,poly> p1,tuple<poly,poly,poly> p2){
            
        }
        tuple<poly,poly,poly> intersect(tuple<poly,poly,poly> p1,tuple<poly,poly,poly> p2){
            return {get<1>(p1) * get<2>(p2) - get<2>(p1) * get<1>(p2), get<2>(p1) * get<0>(p2)-get<0>(p1) * get<2>(p2), get<0>(p1) * get<1>(p2) - get<1>(p1) * get<0>(p2)};
        }
        tuple<poly,poly,poly> parallel_line(tuple<poly,poly,poly> p1,tuple<poly,poly,poly> p2){
            
        }
        tuple<poly,poly,poly,poly,poly,poly> circle_1(tuple<poly,poly,poly> center, double dist){
            
        }
        tuple<poly,poly,poly,poly,poly,poly> circle_2(tuple<poly,poly,poly> p1, tuple<poly,poly,poly> p2, tuple<poly,poly,poly> p3){
            
        }
        void simplify_3(tuple<poly,poly,poly> p){

        }
        tuple<poly,poly,poly> tangent(tuple<poly,poly,poly,poly,poly,poly> circle, tuple<poly,poly,poly> point){

        }
        tuple<poly,poly,poly> area(tuple<poly,poly,poly> p1, tuple<poly,poly,poly> p2, tuple<poly,poly,poly> p3){

        }
};