#include <iostream>

#include "femshell.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    FEMShell windown;
    windown.show();

    return a.exec();

}

