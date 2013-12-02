#ifndef FEMSHELL_H
#define FEMSHELL_H

#include <QMainWindow>

#include "node.h"
#include "elementqn.h"
#include "boundary.h"

namespace Ui {
    class FEMShell;
}

enum FEA{
    ThinPlate,
    ThickPlate,
    Shell
};

enum MeshType{
    Rectangular,
    Circular,
    Ring,
    Cylinder
};



class FEMShell : public QMainWindow
{
        Q_OBJECT

    public:

        Boundary *boundaries;
        ElementQN **elements;
        Node **nodes;

        int nNodes;
        int nElements;


        FEA solver;
        MeshType meshType;

        double E,v,G,K,t,lx,ly,re,ri;
        int nx, ny, nr, npx, npy;



        explicit FEMShell(QWidget *parent = 0);
        ~FEMShell();

        void generateMesh(void);
        void setupRetangularMesh(void);

        public slots:

        void updateSelectedSolverOption(void);
        void updateSelectedMeshOption(void);
        void updateValidOptions(void);
        void updateData(void);

        void updateMeshTypeValidOptions1(void);
        void updateMeshTypeValidOptions2(void);
        void updateMeshTypeValidOptions3(void);
        void updateMeshTypeValidOptions4(void);





    private:
        Ui::FEMShell *ui;
};

#endif // FEMSHELL_H