TEMPLATE = app
CONFIG += console
#CONFIG -= app_bundle
#CONFIG -= qt

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS *= -fopenmp

INCLUDEPATH += ../library/lapacke

LIBS += -L$$PWD/../library/lapacke/ -llapacke -lptlapack -lptcblas -lptf77blas -latlas -lgfortran

#LIBS += -L$$PWD/../library/lapacke/ -llapacke -llapack -lcblas -lf77blas -latlas -lgfortran

LIBS += -lquadmath

LIBS += -lGLU

#INCLUDEPATH += ../library/gnuplot


SOURCES += main.cpp \
    polynomial2d.cpp \
    polynomial1d.cpp \
    node.cpp \
    functor2d.cpp \
    elementdkt.cpp \
    mainwindow.cpp \
    graphicwindow.cpp \
    elementqn.cpp \
    lagrange.cpp \
    integralgauss2d.cpp \
    arcball.cpp \
    elementsqn.cpp \
    femshell.cpp \
    boundary.cpp \
    thickplatemesh.cpp \
    thinplatemesh.cpp \
    element.cpp \
    matrix.cpp \
    flatshellmesh.cpp \
    mesh.cpp \
    elementssqn.cpp \
    shellmesh.cpp \
    gnuplot.cpp \
    thinshellmesh.cpp \
    elementsdkt.cpp

HEADERS += \
    polynomial2d.h \
    polynomial1d.h \
    node.h \
    functor2d.h \
    elementdkt.h \
    mainwindow.h \
    graphicwindow.h \
    mesh.h \
    matrix.h \
    elementqn.h \
    lagrange.h \
    integralgauss2d.h \
    arcball.h \
    elementsqn.h \
    femshell.h \
    boundary.h \
    thickplatemesh.h \
    thinplatemesh.h \
    element.h \
    flatshellmesh.h \
    elementssqn.h \
    shellmesh.h \
    gnuplot.h \
    thinshellmesh.h \
    elementsdkt.h

FORMS += \
    mainwindow.ui \
    femshell.ui



