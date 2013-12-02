#ifndef GRAPHICWINDOW_H
#define GRAPHICWINDOW_H

#include <QGLWidget>
#include <QMouseEvent>
#include <QCursor>

#include "mesh.h"
#include "mesh2.h"
#include "mesh3.h"

#include "arcball.h"


class GraphicWindow : public QGLWidget
{
        Q_OBJECT
    public:
        explicit GraphicWindow(QWidget *parent = 0);

        double xmax, xmin, ymax, ymin, panX, panY;
        bool isMousePress;

        Mesh3 *mesh;

        QVector3D mousePos; //arcball
        QVector2D mousePt; //arcball
        QPoint lastPos; //arcball

        double* transformM;
        double screenW, screenH;

        ArcBall *arcb;

        GLdouble model_view[16];
        GLdouble projection[16];
        GLint viewport[4];

        double phiRot, thetaRot, psiRot;

        QVector3D getMousePos (int x, int y);
        QVector2D normalizeMouse(QPoint qp);

        void euler2matr ();
        void quat2matr (QQuaternion q);
        QVector3D quat2euler (QQuaternion q);
        QQuaternion quatfromEuler ();

    signals:
        
    public slots:

        void initializeGL();
        void resizeGL(int width, int height);
        void paintGL();

        void wheelEvent(QWheelEvent *event);
        void mousePressEvent(QMouseEvent *event);
        void mouseMoveEvent(QMouseEvent *event);
        void mouseReleaseEvent(QMouseEvent *event);
        void mouseDoubleClickEvent(QMouseEvent *event);
        
};

#endif // GRAPHICWINDOW_H


