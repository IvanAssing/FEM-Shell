#ifndef FEMSHELL_H
#define FEMSHELL_H

#include <QMainWindow>


#include "mesh.h"
#include "node.h"
#include "element.h"
#include "boundary.h"

#include "elementdkt.h"
#include "elementsdkt.h"
#include "elementqn.h"
#include "elementsqn.h"
#include "elementssqn.h"

namespace Ui {
    class FEMShell;
}

enum FEA{
    ThinPlate,
    ThickPlate,
    ThinShell,
    ThickShell,
    Solid3D
};

enum MeshType{
    Rectangular,
    Curved,
    Ring,
    Cylinder
};



class FEMShell : public QMainWindow
{
        Q_OBJECT

    public:
        Boundary *boundaries;

        ElementDKT **elementsdkt;
        ElementSDKT **elementssdkt;
        ElementQN **elementsqn;
        ElementSQN **elementssqn;
        ElementSSQN **elementsssqn;

        Node **nodes;
        Mesh *mesh;

        int nNodes;
        int nElements;

        FEA solver;
        MeshType meshType;

        double E,v,G,K,t,lx,ly,re,ri,alpha;
        int nx, ny, npx, npy;

        bool selectiveIntegration;

        explicit FEMShell(QWidget *parent = 0);
        ~FEMShell();

        void readBoundary(void);
        void setupRetangularMesh(void);
        void setupTriangularMesh(void);
        void setupTriangularShellMesh(void);
        void setupCurvedMesh(void);
        void setupRingMesh(void);
        void setupCylinderMesh(void);
        void setupCurvedShellMesh(void);
        void setupRetangularShellMesh(void);
        void setupRingShellMesh(void);



        public slots:

        void updateData(void);
        void updateSelectedMeshOption(void);
        void updateSelectedSolverOption(void);
        void updateGraphic(void);
        void showPlot(void);

        void readTable(int i, int j);
        void updateTable();

        void createMesh(void);
        void solve(void);

        void open(void);
        void save(void);


    private:
        Ui::FEMShell *ui;
};

#endif // FEMSHELL_H
