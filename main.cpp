/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: acolomitchi
 *
 * Created on 28 July 2018, 3:49 PM
 */

#include <cstdlib>

#include <type_traits>
#include <map>
#include <iostream>
#include <string>
#include "geom2d.hpp"


/*
template <class T, typename U, typename R=void> class resolver {
public:
    typedef R (T::*method_t)(const U&) const;
};

template <auto x> struct god_help_me {};

template <class Decl, typename Ret, typename U, Ret (Decl::*fn)(U) const>
struct god_help_me<fn>
{
    god_help_me() {}
    Ret operator()(const Decl& o, U args) {
        return (o.*fn)(args);
    }
};


struct A {
    int incr(const int& x) const { return x+1; }
};
*/

struct sorta_line {
    double sx, sy;
    double ex, ey;
    geomalgos2d::linear_variety type;
    
public:
    sorta_line(
        double x0, double y0, 
        double x1, double y1, 
        geomalgos2d::linear_variety t
    ) :sx(x0), sy(y0), ex(x1), ey(y1),type(t)
    {}
    
    bool intersection_params(const sorta_line& o, double& thisP, double& oP) const
    {
        using linegeom=geomalgos2d::linegeo2d<double, double>;
        return linegeom::intersection_params(
            sx, sy, ex, ey,
            o.sx, o.sy, o.ex, o.ey,
            thisP, oP,
            type, o.type
        );
    }
};

std::ostream& operator<<(std::ostream& o, const sorta_line& l) {
    o << "[" << l.sx <<","<<l.sy<< "-" << l.type << "->"
      << l.ex << "," << l.ey << "]";
    return o;
}

void printIntersection(
    const std::string& what,
    const sorta_line& l0, const sorta_line& l1
) {
    double p0, p1;
    bool has=l0.intersection_params(l1, p0, p1);
    std::cout << what << " " << l0 << " intersects " << l1 << "? "
              << std::boolalpha << has
    ;
    if(has) {
        std::cout << " @{" << p0 << "," <<p1<<"}";
        double x=(1-p0)*l0.sx+p0*l0.ex;
        double y=(1-p0)*l0.sy+p0*l0.ey;
        std::cout<< "["<<x << ","<<y<<"]";
        x=(1-p1)*l1.sx+p1*l1.ex;
        y=(1-p1)*l1.sy+p1*l1.ey;
        std::cout<< "["<<x << ","<<y<<"]";
    }
    std::cout << std::endl;
}
int main(int argc, char** argv) {
//    A a;
//    resolver<A,int,int>::method_t x=&A::incr;
//    god_help_me<&A::incr> y;
//    std::cout << y(a, 2) << std::endl;
//    std::cout << god_help_me<&A::incr>()(a, 2) << std::endl; 
    // -------------------------------^^---LOL
    using linegeom=geomalgos2d::linegeo2d<double, double>;
    using lv=geomalgos2d::linear_variety;
    
    printIntersection(
      "line X line",
      sorta_line(0, 0, 0.25, 0.25, lv::line),
      sorta_line(0, 1, 1, 0, lv::line)
    );
    
    
    // degenerated segments X axis
    printIntersection(
      "degenerated x s-s contained",
      sorta_line(-1, 0, 1, 0, lv::segment),
      sorta_line(-6, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained",
      sorta_line(1, 0, -1, 0, lv::segment),
      sorta_line(-6, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained",
      sorta_line(-1, 0, 1, 0, lv::segment),
      sorta_line(2, 0, -6, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained",
      sorta_line(1, 0, -1, 0, lv::segment),
      sorta_line(2, 0, -6, 0, lv::segment)
    );
    
    printIntersection(
      "degenerated x s-s contained w limit",
      sorta_line(-2, 0, 2, 0, lv::segment),
      sorta_line(-6, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained w limit",
      sorta_line(2, 0, -2, 0, lv::segment),
      sorta_line(-6, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained w limit",
      sorta_line(-2, 0, 2, 0, lv::segment),
      sorta_line(2, 0, -6, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s contained w limit",
      sorta_line(2, 0, -2, 0, lv::segment),
      sorta_line(2, 0, -6, 0, lv::segment)
    );

    // y axis
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, 3, 0, -1, lv::segment),
      sorta_line(0, -3, 0, 1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, -1, 0, 3, lv::segment),
      sorta_line(0, -3, 0, 1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, 3, 0, -1, lv::segment),
      sorta_line(0, 1, 0, -3, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, -1, 0, 3, lv::segment),
      sorta_line(0, 1, 0, -3, lv::segment)
    );
    
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, -3, 0, 1, lv::segment),
      sorta_line(0, 3, 0, -1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, -3, 0, 1, lv::segment),
      sorta_line(0, -1, 0, 3, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, 1, 0, -3, lv::segment),
      sorta_line(0, 3, 0, -1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s overlapping",
      sorta_line(0, 1, 0, -3, lv::segment),
      sorta_line(0, -1, 0, 3, lv::segment)
    );
    
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, -6, 0, 2, lv::segment),
      sorta_line(0, -1, 0, 1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, 2, 0, -6, lv::segment),
      sorta_line(0, -1, 0, 1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, -6, 0, 2, lv::segment),
      sorta_line(0, 1, 0, -1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, 2, 0, -6, lv::segment),
      sorta_line(0, 1, 0, -1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, -1, 0, 1, lv::segment),
      sorta_line(0, -6, 0, 2, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, -1, 0, 1, lv::segment),
      sorta_line(0, 2, 0, -6, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, 1, 0, -1, lv::segment),
      sorta_line(0, -6, 0, 2, lv::segment)
    );
    printIntersection(
      "degenerated y s-s contained",
      sorta_line(0, 1, 0, -1, lv::segment),
      sorta_line(0, 2, 0, -6, lv::segment)
    );
    // degenerated s-s x 1 point
    printIntersection(
      "degenerated x s-s one end",
      sorta_line(0, 0, 1, 0, lv::segment),
      sorta_line(1, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s one end",
      sorta_line(1, 0, 0, 0, lv::segment),
      sorta_line(1, 0, 2, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s one end",
      sorta_line(0, 0, 1, 0, lv::segment),
      sorta_line(2, 0, 1, 0, lv::segment)
    );
    printIntersection(
      "degenerated x s-s one end",
      sorta_line(1, 0, 0, 0, lv::segment),
      sorta_line(2, 0, 1, 0, lv::segment)
    );
    // degenerated s-s y 1 point
    printIntersection(
      "degenerated y s-s one end",
      sorta_line(0, 0, 0, 1, lv::segment),
      sorta_line(0, 1, 0, 2, lv::segment)
    );
    printIntersection(
      "degenerated y s-s one end",
      sorta_line(0, 1, 0, 0, lv::segment),
      sorta_line(0, 1, 0, 2, lv::segment)
    );
    printIntersection(
      "degenerated y s-s one end",
      sorta_line(0, 0, 0, 1, lv::segment),
      sorta_line(0, 2, 0, 1, lv::segment)
    );
    printIntersection(
      "degenerated y s-s one end",
      sorta_line(0, 1, 0, 0, lv::segment),
      sorta_line(0, 2, 0, 1, lv::segment)
    );
    
   
    // ray segment x
    printIntersection(
       "degenerated x r-s",
       sorta_line(0, 0, 1, 0, lv::ray),
       sorta_line(0.5, 0, 0.75, 0, lv::segment)
    );
    printIntersection(
       "degenerated x r-s",
       sorta_line(1, 0, 0, 0, lv::ray),
       sorta_line(0.5, 0, 0.75, 0, lv::segment)
    );
    printIntersection(
       "degenerated x r-s",
       sorta_line(0, 0, 1, 0, lv::ray),
       sorta_line(0.75, 0, 0.5, 0, lv::segment)
    );
    printIntersection(
       "degenerated x r-s",
       sorta_line(1, 0, 0, 0, lv::ray),
       sorta_line(0.75, 0, 0.5, 0, lv::segment)
    );
    
    printIntersection(
       "degenerated x s-r",
       sorta_line(0.5, 0, 0.75, 0, lv::segment),
       sorta_line(0, 0, 1, 0, lv::ray)
    );
    printIntersection(
       "degenerated x s-r",
       sorta_line(0.5, 0, 0.75, 0, lv::segment),
       sorta_line(1, 0, 0, 0, lv::ray)
    );
    printIntersection(
       "degenerated x s-r",
       sorta_line(0.75, 0, 0.5, 0, lv::segment),
       sorta_line(0, 0, 1, 0, lv::ray)
    );
    printIntersection(
       "degenerated x s-r",
       sorta_line(0.75, 0, 0.5, 0, lv::segment),
       sorta_line(1, 0, 0, 0, lv::ray)
    );

    // ray segment y
    printIntersection(
       "degenerated y r-s",
       sorta_line(0, 0, 0, 1, lv::ray),
       sorta_line(0, 0.5, 0, 0.75, lv::segment)
    );
    printIntersection(
       "degenerated y r-s",
       sorta_line(0, 1, 0, 0, lv::ray),
       sorta_line(0, 0.5, 0, 0.75, lv::segment)
    );
    printIntersection(
       "degenerated y r-s",
       sorta_line(0, 0, 0, 1, lv::ray),
       sorta_line(0, 0.75, 0, 0.5, lv::segment)
    );
    printIntersection(
       "degenerated y r-s",
       sorta_line(0, 1, 0, 0, lv::ray),
       sorta_line(0, 0.75, 0, 0.5, lv::segment)
    );
    
    printIntersection(
       "degenerated y s-r",
       sorta_line(0, 0.5, 0, 0.75, lv::segment),
       sorta_line(0, 0, 0, 1, lv::ray)
    );
    printIntersection(
       "degenerated y s-r",
       sorta_line(0, 0.5, 0, 0.75, lv::segment),
       sorta_line(0, 1, 0, 0, lv::ray)
    );
    printIntersection(
       "degenerated y s-r",
       sorta_line(0, 0.75, 0, 0.5, lv::segment),
       sorta_line(0, 0, 0, 1, lv::ray)
    );
    printIntersection(
       "degenerated y s-r",
       sorta_line(0, 0.75, 0, 0.5, lv::segment),
       sorta_line(0, 1, 0, 0, lv::ray)
    );
    
    // ray seg degenerated limits x
    printIntersection(
      "degenerated x s-r one end",
      sorta_line(0, 0, 1, 0, lv::segment),
      sorta_line(0.25, 0, 0.5, 0, lv::ray)
    );
    printIntersection(
      "degenerated x s-r one end",
      sorta_line(1, 0, 0, 0, lv::segment),
      sorta_line(0.25, 0, 0.5, 0, lv::ray)
    );
    printIntersection(
      "degenerated x s-r one end",
      sorta_line(0, 0, 1, 0, lv::segment),
      sorta_line(0.5, 0, 0.25, 0, lv::ray)
    );
    printIntersection(
      "degenerated x s-r one end",
      sorta_line(1, 0, 0, 0, lv::segment),
      sorta_line(0.5, 0, 0.25, 0, lv::ray)
    );
    
    // ray seg degenerated limits y
    printIntersection(
      "degenerated y s-r one end",
      sorta_line(0, 0, 0, 1, lv::segment),
      sorta_line(0, 0.25, 0, 0.5, lv::ray)
    );
    printIntersection(
      "degenerated y s-r one end",
      sorta_line(0, 1, 0, 0, lv::segment),
      sorta_line(0, 0.25, 0, 0.5, lv::ray)
    );
    printIntersection(
      "degenerated y s-r one end",
      sorta_line(0, 0, 0, 1, lv::segment),
      sorta_line(0, 0.5, 0, 0.25, lv::ray)
    );
    printIntersection(
      "degenerated y s-r one end",
      sorta_line(0, 1, 0, 0, lv::segment),
      sorta_line(0, 0.5, 0, 0.25, lv::ray)
    );

    return 0;
}

