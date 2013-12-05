#include "femshell.h"
#include "ui_femshell.h"



#include "thickplatemesh.h"
#include "thinplatemesh.h"


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


    QObject::connect(ui->radioButton, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_2, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_3, SIGNAL(released()), this, SLOT(updateSelectedSolverOption()));
    QObject::connect(ui->radioButton_4, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_5, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_6, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));
    QObject::connect(ui->radioButton_7, SIGNAL(released()), this, SLOT(updateSelectedMeshOption()));

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

    QObject::connect(ui->actionSolver, SIGNAL(triggered()), this, SLOT(solve()));

    //QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

    ui->radioButton_8->setDisabled(true);
    ui->table1->verticalHeader()->setVisible(false);

    boundaries = new Boundary[4];


}


void FEMShell::updateSelectedSolverOption(void)
{
    if(ui->radioButton->isChecked())
    {
        this->solver = ThinPlate;
        ui->radioButton_4->setChecked(true);
        ui->radioButton_5->setDisabled(true);
        ui->radioButton_6->setDisabled(true);
        ui->radioButton_7->setDisabled(true);

        this->meshType = Rectangular;

        ui->lineEdit_13->setDisabled(true); // G
        ui->lineEdit_14->setDisabled(true); // K
        ui->lineEdit_6->setDisabled(false); // ly
        ui->lineEdit_15->setDisabled(true); // re
        ui->lineEdit_16->setDisabled(true); // ri
        ui->lineEdit_9->setDisabled(false); // ny
        ui->lineEdit_18->setDisabled(true); // nr
        ui->lineEdit_10->setDisabled(true); // npx
        ui->lineEdit_17->setDisabled(true); // npy


    }

    if(ui->radioButton_2->isChecked())
    {
        this->solver = ThickPlate;
        ui->radioButton_4->setDisabled(false);
        ui->radioButton_5->setDisabled(false);
        ui->radioButton_6->setDisabled(false);
        ui->radioButton_7->setDisabled(true);
        ui->radioButton_4->setChecked(true);
    }

    if(ui->radioButton_3->isChecked())
    {
        this->solver = Shell;
        ui->radioButton_4->setDisabled(false);
        ui->radioButton_5->setDisabled(false);
        ui->radioButton_6->setDisabled(false);
        ui->radioButton_7->setDisabled(false);
        ui->radioButton_4->setChecked(true);
    }


}

void FEMShell::updateSelectedMeshOption(void)
{

    if(ui->radioButton_4->isChecked())
    {
        this->meshType = Rectangular;
        this->updateMeshTypeValidOptions1();
    }
    if(ui->radioButton_5->isChecked())
    {
        this->meshType = Circular;
        this->updateMeshTypeValidOptions2();
    }
    if(ui->radioButton_6->isChecked())
    {
        this->meshType = Ring;
        this->updateMeshTypeValidOptions3();
    }
    if(ui->radioButton_7->isChecked())
    {
        this->meshType = Cylinder;
        this->updateMeshTypeValidOptions4();
    }


}


void FEMShell::updateValidOptions(void)
{

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
    nx = ui->lineEdit_8->text().toInt();
    ny = ui->lineEdit_9->text().toInt();
    nr = ui->lineEdit_18->text().toInt();
    npx = ui->lineEdit_10->text().toInt();
    npy = ui->lineEdit_17->text().toInt();

}


void FEMShell::updateMeshTypeValidOptions1(void) // Rectangular
{
    ui->lineEdit_13->setDisabled(false); // G
    ui->lineEdit_14->setDisabled(false); // K
    ui->lineEdit_6->setDisabled(false); // ly
    ui->lineEdit_15->setDisabled(!false); // re
    ui->lineEdit_16->setDisabled(!false); // ri
    ui->lineEdit_9->setDisabled(false); // ny
    ui->lineEdit_18->setDisabled(!false); // nr
    ui->lineEdit_10->setDisabled(false); // npx
    ui->lineEdit_17->setDisabled(false); // npy

    this->generateMesh();

}

void FEMShell::updateMeshTypeValidOptions2(void) // Circular
{

    ui->lineEdit_13->setDisabled(false); // G
    ui->lineEdit_14->setDisabled(false); // K
    ui->lineEdit_6->setDisabled(!false); // ly
    ui->lineEdit_15->setDisabled(false); // re
    ui->lineEdit_16->setDisabled(false); // ri
    ui->lineEdit_9->setDisabled(!false); // ny
    ui->lineEdit_18->setDisabled(!false); // nr
    ui->lineEdit_10->setDisabled(false); // npx
    ui->lineEdit_17->setDisabled(false); // npy

}

void FEMShell::updateMeshTypeValidOptions3(void) // Ring
{

    ui->lineEdit_13->setDisabled(false); // G
    ui->lineEdit_14->setDisabled(false); // K
    ui->lineEdit_6->setDisabled(!false); // ly
    ui->lineEdit_15->setDisabled(false); // re
    ui->lineEdit_16->setDisabled(false); // ri
    ui->lineEdit_9->setDisabled(!false); // ny
    ui->lineEdit_18->setDisabled(false); // nr
    ui->lineEdit_10->setDisabled(false); // npx
    ui->lineEdit_17->setDisabled(false); // npy

}


void FEMShell::updateMeshTypeValidOptions4(void) // Cylinder
{

    ui->lineEdit_13->setDisabled(false); // G
    ui->lineEdit_14->setDisabled(false); // K
    ui->lineEdit_6->setDisabled(!false); // ly
    ui->lineEdit_15->setDisabled(false); // re
    ui->lineEdit_16->setDisabled(false); // ri
    ui->lineEdit_9->setDisabled(!false); // ny
    ui->lineEdit_18->setDisabled(false); // nr
    ui->lineEdit_10->setDisabled(false); // npx
    ui->lineEdit_17->setDisabled(false); // npy

}

void FEMShell::solve(void)
{
    this->updateTable();
    mesh->solve();
}


FEMShell::~FEMShell()
{
    delete ui;
}

void FEMShell::generateMesh(void)
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

    nx = ny = 20;
    npx = npy = 3;
    lx = ly = 1;



    double vi = 0.3;
    double E = 200e9;
    double t = 0.02;
    double G = 75e9;

    double GKt = G*5./6.*t;

    double Ept = E*t*t*t/(12.*(1.0-vi*vi));
    double Em = E/(1.0-vi*vi);

    Matrix D(3,3);

    D(0, 0) = Ept;
    D(0, 1) = Ept*vi;
    D(1, 0) = Ept*vi;
    D(1, 1) = Ept;
    D(2, 2) = Ept*(1-vi)/2.0;

    Matrix Dm(3,3);

    Dm(0, 0) = Em;
    Dm(0, 1) = Em*vi;
    Dm(1, 0) = Em*vi;
    Dm(1, 1) = Em;
    Dm(2, 2) = Em*(1-vi)/2.0;


    //this->setupRetangularMesh();
    //this->setupTriangularMesh();

    ri = 4.0;
    re = 5.0;
    alpha = M_PI;

    //this->setupCurvedMesh();
    this->setupRingMesh();

    this->mesh = new ThickPlateMesh(nNodes, nodes, nElements, elementsqn, npx-1, npy-1, D, GKt);

    //this->mesh = new FlatShel (nNodes, nodes, nElements, elementsqn, 1, 1, D, GKt, Dm);

    //this->mesh = new ThinPlateMesh(nNodes, nodes, nElements, elementsdkt, D);


    ui->widget->mesh = this->mesh;

}


void FEMShell::setupRetangularMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    int nNodesElement = npx*npy;

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

    Node **ptrNodes = new Node*[nNodesElement];

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
            elementsqn[elementIndex++] = new ElementQN(nNodesElement, ptrNodes);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
}

void FEMShell::setupRingMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1);

    int nNodesElement = npx*npy;

    nNodes = nNodesx*nNodesy;

    nodes = new Node*[nNodes];

    double dx = (re-ri)/(nNodesx-1);
    double dy = 2.0*M_PI/(nNodesy);


    for(int i=0; i<nNodesy; i++)
        for(int j=0; j<nNodesx; j++)
            nodes[j+nNodesx*i] = new Node(j+nNodesx*i, (ri+j*dx)*cos(dy*i), (ri+j*dx)*sin(dy*i));


//    for(int i=0; i<nNodesx; i++)
//        nodes[i]->setup(boundaries[0]);

//    for(int i=nNodes-nNodesx; i<nNodes; i++)
//        nodes[i]->setup(boundaries[1]);

        for(int i=0; i<nNodesy; i++)
            nodes[i*nNodesx]->setup(boundaries[2]);

        for(int i=1; i<=nNodesy; i++)
            nodes[i*nNodesx-1]->setup(boundaries[3]);

    Node **ptrNodes = new Node*[nNodesElement];

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
            elementsqn[elementIndex++] = new ElementQN(nNodesElement, ptrNodes);
        }

    for(int je=0; je<nx; je++)
    {
        ni = 0;
        for(int i=0; i<npy-1; i++)
            for(int j=0; j<npx; j++)
                ptrNodes[ni++] = nodes[(ny-1)*(npx-1)*nNodesx + i*nNodesx + je*(npx-1) + j];
        for(int j=0; j<npx; j++)
            ptrNodes[ni++] = nodes[je*(npx-1) + j];
        elementsqn[elementIndex++] = new ElementQN(nNodesElement, ptrNodes);
    }



    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));

}

void FEMShell::setupCurvedMesh(void)
{

    int nNodesx = nx*(npx-1)+1;
    int nNodesy = ny*(npy-1)+1;

    int nNodesElement = npx*npy;

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

    Node **ptrNodes = new Node*[nNodesElement];

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
            elementsqn[elementIndex++] = new ElementQN(nNodesElement, ptrNodes);
        }

    ui->table1->setRowCount(nNodes);

    this->updateTable();

    QObject::connect(ui->table1, SIGNAL(cellChanged(int,int)), this, SLOT(readTable(int,int)));
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

void FEMShell::setupBoundaryConditions(void)
{

}

