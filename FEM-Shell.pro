TEMPLATE = app
CONFIG += console
#CONFIG -= app_bundle
#CONFIG -= qt

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS *= -fopenmp

INCLUDEPATH += ../library/lapack

LIBS += -L$$PWD/../library/lapack/ -llapack_linux_x86-64 -lf2c -lcblaswr -lptcblas -latlas

LIBS += -lquadmath

LIBS += -lGLU


SOURCES += main.cpp \
    polynomial2d.cpp \
    polynomial1d.cpp \
    node.cpp \
    functor2d.cpp \
    elementdkt.cpp \
    mainwindow.cpp \
    graphicwindow.cpp \
    mesh.cpp \
    matrix.cpp \
    elementqn.cpp \
    mesh2.cpp \
    lagrange.cpp \
    integralgauss2d.cpp \
    arcball.cpp

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
    mesh2.h \
    lagrange.h \
    integralgauss2d.h \
    arcball.h

FORMS += \
    mainwindow.ui



