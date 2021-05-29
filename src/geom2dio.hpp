#ifndef GEOM2DIO_HPP
#define GEOM2DIO_HPP

#include <iostream>
#include <map>
#include <string>

#include "linear2d.hpp"

namespace geomalgos2d
{
std::ostream& operator<<(
    std::ostream& out,
    const linear_variety value
) {
    static std::map<linear_variety, std::string> strings;
    if (strings.size() == 0){
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(ray);
        INSERT_ELEMENT(segment);
        INSERT_ELEMENT(line);
#undef INSERT_ELEMENT
    }
    return out << strings[value];
}

} // geomalgos2d
#endif // GEOM2DIO_HPP
