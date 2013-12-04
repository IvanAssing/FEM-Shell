#include <iostream>

#include "node.cpp"
#include "elementdkt.h"


#include "mainwindow.h"
#include "femshell.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);


    FEMShell windown;

    windown.show();

    return a.exec();

}

