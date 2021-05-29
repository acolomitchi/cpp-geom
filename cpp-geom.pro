TEMPLATE = app
TARGET = geom2d-tests
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    tests/main.cpp \
    tests/line-intersect.cpp \
    tests/closest-point.cpp

HEADERS += \
    src/prec2d.hpp \
    src/linear2d.hpp \
    src/geom2dio.hpp \
    src/geom2d.hpp \
    src/dist2d.hpp \
    tests/sorta_line.hpp

INCLUDEPATH += ../../../Catch2/single_include
