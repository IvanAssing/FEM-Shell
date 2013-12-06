#include "femshell.h"
#include "ui_femshell.h"

#include <QMessageBox>

#include "thickplatemesh.h"
#include "thinplatemesh.h"
#include "flatshellmesh.h"
#include "shellmesh.h"


FEMShell::FEMShell(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::FEMShell)
{
    ui->setupUi(this);


    ui->table1->setColumnCount(16);
    ui->table1->setHorizontalHeaderLabels(
                QStringList() << "index" << "x" << "y" << "z" <<"Lock x" << "Fx or u"
                <<"Lock y" << "Fy or v" <<"Lock z" << "Fz or w" <<"Lock θx" << "Mx or θx"
                <<"Lock θy" << "My or θy" <<"Lock θz" << "Mz or θz"
                );

    QObject::connect(ui->lineEdit_5, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_6, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_7, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_8, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_9, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_10, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_11, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_12, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_13, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_14, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_15, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_16, SIGNAL(editingFinished()), this, SLOT(updateData()));
    QObject::connect(ui->lineEdit_17, SIGNAL(editingFinished()), this, SLOT(updateData()));

    QObject::connect(ui->radioButton, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_2, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_3, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_4, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_5, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_6, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_7, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));

    QObject::connect(ui->radioButton_9, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_10, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_11, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_12, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_13, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_14, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_15, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_16, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_17, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_18, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->radioButton_19, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->checkBox_26, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->checkBox_27, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->checkBox_28, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->checkBox_29, SIGNAL(released()), this, SLOT(updateGraphic()));
    QObject::connect(ui->lineEdit_18, SIGNAL(returnPressed()), this, SLOT(updateGraphic()));

    QObject::connect(ui->actionSolver, SIGNAL(triggered()), this, SLOT(solve()));

    QObject::connect(ui->action_Genarate_Mesh, SIGNAL(triggered()), this, SLOT(createMesh()));

    ui->radioButton_8->setDisabled(true);
    ui->table1->verticalHeader()->setVisible(false);

    boundaries = new Boundary[4];

    mesh = NULL;


}


void FEMShell::createMesh(void)
{
    updateData();
    readBoundary();

//    v = 0.3;
//    E = 200e9;
//    t = 0.02;
//    G = 75e9;
//    K = 5./6.;

//    ny = 10;
//    nx = 10;
//    npx = npy = 2;
//    lx = ly = 1;
//    ri = 4.0;
//    re = 5.0;
//    alpha = M_PI/4;

    double GKt = G*K*t;

    double Ef = E*t*t*t/(12.*(1.0-v*v));
    double Em = E/(1.0-v*v);

    Matrix D(3,3);

    D(0, 0) = Ef;
    D(0, 1) = Ef*v;
    D(1, 0) = Ef*v;
    D(1, 1) = Ef;
    D(2, 2) = Ef*(1-v)/2.0;

    Matrix Dm(3,3);

    Dm(0, 0) = Em;
    Dm(0, 1) = Em*v;
    Dm(1, 0) = Em*v;
    Dm(1, 1) = Em;
    Dm(2, 2) = Em*(1-v)/2.0;

    if(solver == ThinPlate && meshType == Rectangular)
    {
        setupTriangularMesh();
        this->mesh = new ThinPlateMesh(nNodes, nodes, nElements, elementsdkt, D);
        ui->widget->mesh = this->mesh;

    }

    if(solver == ThickPlate)
    {
        if(meshType == Rectangular)
        {
            setupRetangularMesh();
            this->mesh = new ThickPlateMesh(nNodes, nodes, nElements, elementsqn, npx, npy, D, GKt);
            ui->widget->mesh = this->mesh;
        }

        if(meshType == Curved)
        {
            setupCurvedMesh();
            this->mesh = new ThickPlateMesh(nNodes, nodes, nElements, elementsqn, npx, npy, D, GKt);
            ui->widget->mesh = this->mesh;
        }

        if(meshType == Ring)
        {
            setupRingMesh();
            this->mesh = new ThickPlateMesh(nNodes, nodes, nElements, elementsqn, npx, npy, D, GKt);
            ui->widget->mesh = this->mesh;
        }

    }

    if(solver == Shell)
    {
        if(meshType == Rectangular)
        {
            setupRetangularShellMesh();
            this->mesh = new FlatShellMesh(nNodes, nodes, nElements, elementssqn, npx, npy, D, GKt, Dm);
            ui->widget->mesh = this->mesh;
        }

        if(meshType == Curved)
        {
            setupCurvedShellMesh();
            this->mesh = new FlatShellMesh(nNodes, nodes, nElements, elementssqn, npx, npy, D, GKt, Dm);
            ui->widget->mesh = this->mesh;
        }

        if(meshType == Ring)
        {
            setupRingShellMesh();
            this->mesh = new FlatShellMesh(nNodes, nodes, nElements, elementssqn, npx, npy, D, GKt, Dm);
            ui->widget->mesh = this->mesh;
        }

        if(meshType == Cylinder)
        {
            setupCylinderMesh();
            this->mesh = new ShellMesh(nNodes, nodes, nElements, elementsssqn, npx, npy, D, GKt, Dm);
            ui->widget->mesh = this->mesh;
        }

    }


}



void FEMShell::readTable(int i, int j)
{
    QObject::disconnect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

    nodes[i]->x = ui->table1->item(i, 1)->text().toDouble();
    nodes[i]->y = ui->table1->item(i, 2)->text().toDouble();
    nodes[i]->z = ui->table1->item(i, 3)->text().toDouble();

    for(j=0; j<6; j++)
    {
        nodes[i]->lockStatus[j] = ui->table1->item(i, 4+2*j)->text().toInt() == 0 ? false : true;
        nodes[i]->loadValues[j] = ui->table1->item(i, 5+2*j)->text().toDouble();
    }


    ui->table1->setItem(i ,0, new QTableWidgetItem(QString("%1").arg(nodes[i]->index)));
    ui->table1->setItem(i ,1, new QTableWidgetItem(QString("%1").arg(nodes[i]->x, 0, 'e', 3)));
    ui->table1->setItem(i ,2, new QTableWidgetItem(QString("%1").arg(nodes[i]->y, 0, 'e', 3)));
    ui->table1->setItem(i ,3, new QTableWidgetItem(QString("%1").arg(nodes[i]->z, 0, 'e', 3)));
    for(int j=0; j<6; j++)
    {
        ui->table1->setItem(i ,4+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->lockStatus[j])));
        ui->table1->setItem(i ,5+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->loadValues[j], 0, 'e', 3)));
    }

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

}


void FEMShell::updateTable(void)
{
    QObject::disconnect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

    for(int i=0; i<nNodes; i++)
    {
        ui->table1->setItem(i ,0, new QTableWidgetItem(QString("%1").arg(nodes[i]->index)));
        ui->table1->setItem(i ,1, new QTableWidgetItem(QString("%1").arg(nodes[i]->x, 0, 'e', 3)));
        ui->table1->setItem(i ,2, new QTableWidgetItem(QString("%1").arg(nodes[i]->y, 0, 'e', 3)));
        ui->table1->setItem(i ,3, new QTableWidgetItem(QString("%1").arg(nodes[i]->z, 0, 'e', 3)));
        for(int j=0; j<6; j++)
        {
            ui->table1->setItem(i ,4+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->lockStatus[j])));
            ui->table1->setItem(i ,5+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->loadValues[j], 0, 'e', 3)));
        }

    }

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}


void FEMShell::updateData(void)
{

    E = ui->lineEdit_12->text().toDouble();
    v = ui->lineEdit_11->text().toDouble();
    G = ui->lineEdit_13->text().toDouble();
    K = ui->lineEdit_14->text().toDouble();
    t = ui->lineEdit_5->text().toDouble();
    lx = ui->lineEdit_7->text().toDouble();
    ly = ui->lineEdit_6->text().toDouble();
    re = ui->lineEdit_15->text().toDouble();
    ri = ui->lineEdit_16->text().toDouble();
    alpha = M_PI*(ui->lineEdit_16->text().toDouble())/180.;
    nx = ui->lineEdit_8->text().toInt();
    ny = ui->lineEdit_9->text().toInt();
    npx = 1 + ui->lineEdit_10->text().toInt();
    npy = 1 + ui->lineEdit_17->text().toInt();

    selectiveIntegration = ui->checkBox_25->isChecked();

}


void FEMShell::solve(void)
{
    this->updateTable();
    if(mesh)
    {
        int ndof = solver == Shell ? 6*nNodes : 3*nNodes;
        int resp = QMessageBox::warning(this, QString("Solver"),
                                        QString("NDOF = %1 \n Continue?").arg(ndof),
                                        QMessageBox::Ok | QMessageBox::Cancel);

        if(resp == QMessageBox::Ok)
        {
            mesh->solve();
            mesh->plot();

        }



    }
}


FEMShell::~FEMShell()
{
    delete ui;
}


void FEMShell::readBoundary(void)
{
    bool locked[6];
    double load[6];

    locked[0] = ui->checkBox->isChecked();
    locked[1] = ui->checkBox_2->isChecked();
    locked[2] = ui->checkBox_3->isChecked();
    locked[3] = ui->checkBox_4->isChecked();
    locked[4] = ui->checkBox_5->isChecked();
    locked[5] = ui->checkBox_6->isChecked();

    load[0] = ui->lineEdit->text().toDouble();
    load[1] = ui->lineEdit_2->text().toDouble();
    load[2] = ui->lineEdit_3->text().toDouble();
    load[3] = ui->lineEdit_4->text().toDouble();
    load[4] = ui->lineEdit_19->text().toDouble();
    load[5] = ui->lineEdit_20->text().toDouble();

    boundaries[0] = Boundary(0, locked, load);

    locked[0] = ui->checkBox_7->isChecked();
    locked[1] = ui->checkBox_8->isChecked();
    locked[2] = ui->checkBox_9->isChecked();
    locked[3] = ui->checkBox_10->isChecked();
    locked[4] = ui->checkBox_11->isChecked();
    locked[5] = ui->checkBox_12->isChecked();

    load[0] = ui->lineEdit_21->text().toDouble();
    load[1] = ui->lineEdit_22->text().toDouble();
    load[2] = ui->lineEdit_23->text().toDouble();
    load[3] = ui->lineEdit_24->text().toDouble();
    load[4] = ui->lineEdit_25->text().toDouble();
    load[5] = ui->lineEdit_26->text().toDouble();

    boundaries[1] = Boundary(1, locked, load);

    locked[0] = ui->checkBox_13->isChecked();
    locked[1] = ui->checkBox_14->isChecked();
    locked[2] = ui->checkBox_15->isChecked();
    locked[3] = ui->checkBox_16->isChecked();
    locked[4] = ui->checkBox_17->isChecked();
    locked[5] = ui->checkBox_18->isChecked();

    load[0] = ui->lineEdit_27->text().toDouble();
    load[1] = ui->lineEdit_28->text().toDouble();
    load[2] = ui->lineEdit_29->text().toDouble();
    load[3] = ui->lineEdit_30->text().toDouble();
    load[4] = ui->lineEdit_31->text().toDouble();
    load[5] = ui->lineEdit_32->text().toDouble();

    boundaries[2] = Boundary(2, locked, load);

    locked[0] = ui->checkBox_19->isChecked();
    locked[1] = ui->checkBox_20->isChecked();
    locked[2] = ui->checkBox_21->isChecked();
    locked[3] = ui->checkBox_22->isChecked();
    locked[4] = ui->checkBox_23->isChecked();
    locked[5] = ui->checkBox_24->isChecked();

    load[0] = ui->lineEdit_33->text().toDouble();
    load[1] = ui->lineEdit_34->text().toDouble();
    load[2] = ui->lineEdit_35->text().toDouble();
    load[3] = ui->lineEdit_36->text().toDouble();
    load[4] = ui->lineEdit_37->text().toDouble();
    load[5] = ui->lineEdit_38->text().toDouble();

    boundaries[3] = Boundary(3, locked, load);

}


void FEMShell::setupTriangularMesh(void)
{
    int nNodesx = nx+1;
    int nNodesy = ny+1;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = lx/nx;
    double dy = ly/ny;

    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, j*dx, i*dy);

    nElements = 2*nx*ny;
    elementsdkt = new ElementDKT*[nElements];

    int elementIndex = 0;

    for(int i=0; i<ny; i++)
        for(int j=0; j<nx; j++)
        {
            elementsdkt[elementIndex++] = new ElementDKT(elementIndex, nodes[nNodesx*i + j], nodes[nNodesx*i + j + 1], nodes[nNodesx*(i+1) + j + 1]);
            elementsdkt[elementIndex++] = new ElementDKT(elementIndex, nodes[nNodesx*i + j], nodes[nNodesx*(i+1) + j + 1], nodes[nNodesx*(i+1) + j]);

            //            elementsdkt[elementIndex++] = new ElementDKT(elementIndex, nodes[nx*i + j], nodes[nx*i + j + 1], nodes[nx*(i+1) + j]);
            //            elementsdkt[elementIndex++] = new ElementDKT(elementIndex, nodes[nx*i + j + 1], nodes[nx*(i+1) + j + 1], nodes[nx*(i+1) + j]);
        }

    for(int i=0; i<nNodesx; i++)
        nodes[i]->setup(boundaries[0]);

    for(int i=nNodes-nNodesx; i<nNodes; i++)
        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    ui->table1->setRowCount(nNodes);

    for(int i=0; i<nNodes; i++)
    {
        ui->table1->setItem(i ,0, new QTableWidgetItem(QString("%1").arg(nodes[i]->index)));
        ui->table1->setItem(i ,1, new QTableWidgetItem(QString("%1").arg(nodes[i]->x, 0, 'e', 3)));
        ui->table1->setItem(i ,2, new QTableWidgetItem(QString("%1").arg(nodes[i]->y, 0, 'e', 3)));
        ui->table1->setItem(i ,3, new QTableWidgetItem(QString("%1").arg(nodes[i]->z, 0, 'e', 3)));
        for(int j=0; j<6; j++)
        {
            ui->table1->setItem(i ,4+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->lockStatus[j])));
            ui->table1->setItem(i ,5+2*j, new QTableWidgetItem(QString("%1").arg(nodes[i]->loadValues[j], 0, 'e', 3)));
        }

    }

}


void FEMShell::setupRetangularMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = lx/(nNodesx-1);
    double dy = ly/(nNodesy-1);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, j*dx, i*dy);


    for(int i=0; i<nNodesx; i++)
        nodes[i]->setup(boundaries[0]);

    for(int i=nNodes-nNodesx; i<nNodes; i++)
        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementsqn = new ElementQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementsqn[elementIndex++] = new ElementQN(npx, npy, ptrNodes, selectiveIntegration);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}


void FEMShell::setupRingMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1);

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = (re-ri)/(nNodesx-1);
    double dy = 2.0*M_PI/(nNodesy);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, (ri+j*dx)*cos(dy*i), (ri+j*dx)*sin(dy*i));

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementsqn = new ElementQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny-1; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementsqn[elementIndex++] = new ElementQN(npx, npy, ptrNodes, selectiveIntegration);
        }

    for(int je=0; je<nx; je++)
    {
        ni = 0;
        for(int i=0; i<npy-1; i++)
            for(int j=0; j<npx; j++)
                ptrNodes[ni++] = nodes[(ny-1)*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
        for(int j=0; j<npx; j++)
            ptrNodes[ni++] = nodes[je*(npx-1) + j];
        elementsqn[elementIndex++] = new ElementQN(npx, npy, ptrNodes, selectiveIntegration);
    }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

}


void FEMShell::setupCurvedMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = (re-ri)/(nNodesx-1);
    double dy = alpha/(nNodesy-1);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, (ri+j*dx)*cos(-0.5*alpha+dy*i), (ri+j*dx)*sin(-0.5*alpha+dy*i));

    for(int i=0; i<nNodesx; i++)
        nodes[i]->setup(boundaries[0]);

    for(int i=nNodes-nNodesx; i<nNodes; i++)
        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementsqn = new ElementQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementsqn[elementIndex++] = new ElementQN(npx, npy, ptrNodes, selectiveIntegration);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}


void FEMShell::setupCylinderMesh(void)
{
    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1);

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = lx/(nNodesx-1);
    double dy = 2.0*M_PI/(nNodesy);


    for(int i=0; i<nNodesx; i++)
        for(int j=0; j<nNodesy; j++)
            nodes[j+nNodesy*i] = new Node(j+nNodesy*i, i*dx, re*cos(dy*j), re*sin(dy*j));




    //    for(int i=0; i<nNodesx; i++)
    //        nodes[i]->setup(boundaries[0]);

    //    for(int i=nNodes-nNodesx; i<nNodes; i++)
    //        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i]->setup(boundaries[2]);

    for(int i=nNodes-nNodesy; i<nNodes; i++)
        nodes[i]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementsssqn = new ElementSSQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<nx; ie++)
    {
        for(int je=0; je<ny-1; je++)
        {
            ni = 0;
            for(int i=0; i<npx; i++)
                for(int j=0; j<npy; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesy + i*nNodesy + je*(npy-1) + j];
            elementsssqn[elementIndex++] = new ElementSSQN(npx, npy, ptrNodes);
        }

        ni = 0;
        for(int i=0; i<npx; i++)
        {
            for(int j=0; j<npy-1; j++)
                ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesy + i*nNodesy + (ny-1)*(npy-1)+ j];
            ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesy + i*nNodesy];
        }
        elementsssqn[elementIndex++] = new ElementSSQN(npx, npy, ptrNodes);

    }



    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

}


void FEMShell::setupRetangularShellMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = lx/(nNodesx-1);
    double dy = ly/(nNodesy-1);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, j*dx, i*dy);


    for(int i=0; i<nNodesx; i++)
        nodes[i]->setup(boundaries[0]);

    for(int i=nNodes-nNodesx; i<nNodes; i++)
        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementssqn = new ElementSQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementssqn[elementIndex++] = new ElementSQN(npx, npy, ptrNodes, selectiveIntegration);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}


void FEMShell::setupCurvedShellMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = (re-ri)/(nNodesx-1);
    double dy = alpha/(nNodesy-1);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i,
                                          (ri+j*dx)*cos(-0.5*alpha+dy*i), (ri+j*dx)*sin(-0.5*alpha+dy*i));


    for(int i=0; i<nNodesx; i++)
        nodes[i]->setup(boundaries[0]);

    for(int i=nNodes-nNodesx; i<nNodes; i++)
        nodes[i]->setup(boundaries[1]);

    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementssqn = new ElementSQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementssqn[elementIndex++] = new ElementSQN(npx, npy, ptrNodes);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}


void FEMShell::setupRingShellMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1);

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = (re-ri)/(nNodesx-1);
    double dy = 2.0*M_PI/(nNodesy);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, (ri+j*dx)*cos(dy*i), (ri+j*dx)*sin(dy*i));


    for(int i=0; i<nNodesy; i++)
        nodes[i*nNodesx]->setup(boundaries[2]);

    for(int i=1; i<=nNodesy; i++)
        nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[npx*npy];

    nElements = nx*ny;
    elementssqn = new ElementSQN*[nElements];
    int elementIndex = 0;

    int ni = 0;
    for(int ie=0; ie<ny-1; ie++)
        for(int je=0; je<nx; je++)
        {
            ni = 0;
            for(int i=0; i<npy; i++)
                for(int j=0; j<npx; j++)
                    ptrNodes[ni++] = nodes[ie*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
            elementssqn[elementIndex++] = new ElementSQN(npx, npy, ptrNodes, selectiveIntegration);
        }

    for(int je=0; je<nx; je++)
    {
        ni = 0;
        for(int i=0; i<npy-1; i++)
            for(int j=0; j<npx; j++)
                ptrNodes[ni++] = nodes[(ny-1)*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
        for(int j=0; j<npx; j++)
            ptrNodes[ni++] = nodes[je*(npx-1) + j];
        elementssqn[elementIndex++] = new ElementSQN(npx, npy, ptrNodes, selectiveIntegration);
    }



    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

}


void FEMShell::updateSelectedSolverOption(void)
{
    if(ui->radioButton->isChecked())
        this->solver = ThinPlate;

    if(ui->radioButton_2->isChecked())
        this->solver = ThickPlate;

    if(ui->radioButton_3->isChecked())
        this->solver = Shell;

}


void FEMShell::updateSelectedMeshOption(void)
{

    if(ui->radioButton_4->isChecked())
        this->meshType = Rectangular;

    if(ui->radioButton_5->isChecked())
        this->meshType = Curved;

    if(ui->radioButton_6->isChecked())
        this->meshType = Ring;

    if(ui->radioButton_7->isChecked())
        this->meshType = Cylinder;


}


void FEMShell::updateGraphic(void)
{
    if(ui->radioButton_9->isChecked())
        ui->widget->data.var = U;
    if(ui->radioButton_10->isChecked())
        ui->widget->data.var = V;
    if(ui->radioButton_11->isChecked())
        ui->widget->data.var = W;
    if(ui->radioButton_12->isChecked())
        ui->widget->data.var = RX;
    if(ui->radioButton_13->isChecked())
        ui->widget->data.var = RY;
    if(ui->radioButton_14->isChecked())
        ui->widget->data.var = RZ;
    if(ui->radioButton_15->isChecked())
        ui->widget->data.var = MX;
    if(ui->radioButton_16->isChecked())
        ui->widget->data.var = MY;
    if(ui->radioButton_17->isChecked())
        ui->widget->data.var = MXY;
    if(ui->radioButton_18->isChecked())
        ui->widget->data.var = QX;
    if(ui->radioButton_19->isChecked())
        ui->widget->data.var = QY;

    ui->widget->data.nodes = ui->checkBox_26->isChecked();
    ui->widget->data.elements = ui->checkBox_27->isChecked();
    ui->widget->data.load = ui->checkBox_29->isChecked();
    ui->widget->data.def = ui->checkBox_28->isChecked();
    ui->widget->data.factor = ui->lineEdit_18->text().toDouble();

    if(ui->widget->data.factor == 0.0)
        ui->widget->data.factor = 1.0;

    ui->widget->updateGL();
}
