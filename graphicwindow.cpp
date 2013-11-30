#include "graphicwindow.h"

#include <cmath>
#include <iomanip>

#include <QImage>
#include <QString>
#include <QImageWriter>
#include <QDateTime>

#include <GL/glu.h>

static const double deg2rad    =  0.01745329251994329576;
static const double pi_2       =  1.57079632679489661922;
static const double rad2deg    = 57.29577951308232087721;


GraphicWindow::GraphicWindow(QWidget *parent) :
    QGLWidget(parent)
{

    xmin = -10.;
    ymin = -10.;
    xmax = +10.;
    ymax = +10.;

    mesh = new Mesh2;


    this->phiRot           = 0.;
    this->thetaRot         = 0.;
    this->psiRot           = 0.;
    this->mousePos         = QVector3D(0., 0., 0.);

    arcb       = new ArcBall();
    transformM = new double[16];

    for (unsigned int i=1; i<15; ++i)
        transformM[i] = 0.;

    for (unsigned int i=0; i<16; i+=5)
        transformM[i] = 1.;
}


void GraphicWindow::initializeGL()
{
    glShadeModel(GL_SMOOTH);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClearDepth(1.0f);

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

}

void GraphicWindow::resizeGL(int width, int height)
{

    double width_ = static_cast<double>(this->width());
    double height_ = static_cast<double>(this->height());

    if (width_ > height_)
    {
        height_ = height_?height_:1;
        double correction = 0.5 * ( width_/ height_ * (ymax-ymin) - (xmax-xmin) );
        xmin   -= correction;
        xmax +=correction;
    }
    else
    {
        width_ = width_?width_:1;
        double correction = 0.5 * ( height_ / width_ * (xmax-xmin) - (ymax-ymin) );
        ymin   -= correction;
        ymax  += correction;
    }

    glViewport( 0, 0, (GLint)width, (GLint)height );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(xmin, xmax, ymin, ymax, 10.0, -10.0);
    //gluPerspective(60, (float)width/height, 0.1, 50000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    screenW = double(width);
    screenH = double(height);

    arcb->set_bounds (width, height);

}


void GraphicWindow::wheelEvent(QWheelEvent *event)
{
    double zoom_factor = -event->delta()/120*0.2;


    double X_opengl = event->x()/static_cast<double>(this->width())*(xmax - xmin)+xmin;
    double Y_opengl  = (this->height()-event->y())/static_cast<double>(this->height())*(ymax - ymin)+ymin;

    xmin -= (X_opengl-xmin)*zoom_factor;
    xmax += (xmax-X_opengl)*zoom_factor;

    ymin -= (Y_opengl-ymin)*zoom_factor;
    ymax += (ymax-Y_opengl)*zoom_factor;


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(xmin, xmax, ymin, ymax, 10.0, -10.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    this->repaint();

}

void GraphicWindow::mousePressEvent(QMouseEvent *event)
{
    isMousePress = true;
    panX = event->x();
    panY = event->y();

    lastPos = event->pos();
    mousePt = normalizeMouse(lastPos);


    if (event->buttons() & Qt::LeftButton)
    {
        mousePos = getMousePos(event->x(), event->y());

        // Update start vector and prepare
        // For dragging
        arcb->click (mousePt);
    }


}


void GraphicWindow::mouseReleaseEvent(QMouseEvent *event)
{
    isMousePress = false;
}

void GraphicWindow::mouseDoubleClickEvent(QMouseEvent *event)
{
    QDateTime now = QDateTime::currentDateTime();

    QString filename = QString("pics/FEM-Shell-snapshot-")
            + now.toString("yyyyMMddhhmmsszzz") + QString(".png");

    this->updateGL();

    this->grabFrameBuffer(true).save(filename, "PNG", 100);

}


void GraphicWindow::mouseMoveEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::LeftButton)
    {

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        double X_opengl = (-event->x()+panX)/static_cast<double>(this->width())*(xmax - xmin);
        double Y_opengl  = (event->y()-panY)/static_cast<double>(this->height())*(ymax - ymin);

        xmax += X_opengl;
        xmin += X_opengl;

        ymax += Y_opengl;
        ymin += Y_opengl;

        panX = event->x();
        panY = event->y();


        glOrtho(xmin, xmax, ymin, ymax, 1.0, -1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        this->repaint();

    }

    lastPos = event->pos();
    mousePt = normalizeMouse(lastPos);

    if (event->buttons() & Qt::RightButton)
    {
        arcb->drag(mousePt);

        quat2matr(arcb->quaternion);
        QVector3D result = quat2euler(arcb->quaternion)*rad2deg;

        if ( thetaRot!=result.x() )
        {
            thetaRot = result.x();
            if ( thetaRot<-180. ) thetaRot += 360.;
            if ( thetaRot>180. ) thetaRot -= 360.;
            //emit SIG_thetaRotationChanged(int(thetaRot));
        }
        if ( phiRot!=result.z() )
        {
            phiRot = result.z();
            if ( phiRot<-90. ) phiRot += 180.;
            if ( phiRot>90. ) phiRot -= 180.;
            //emit SIG_phiRotationChanged(int(phiRot));
        }
        if ( psiRot!=result.y() )
        {
            psiRot = result.y();
            if ( psiRot<-180. ) psiRot += 180.;
            if ( psiRot>180. ) psiRot -= 180.;
            //emit SIG_psiRotationChanged(int(psiRot));
        }




        updateGL();

    }
}



void GraphicWindow::paintGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glGetDoublev (GL_MODELVIEW_MATRIX, model_view);
    glGetDoublev (GL_PROJECTION_MATRIX, projection);
    glGetIntegerv (GL_VIEWPORT, viewport);

    glMultMatrixd (transformM);



    for(int i=0; i<mesh->nElements; i++)
        mesh->elements[i]->draw();


    for(int i=0; i<mesh->nNodes; i++)
        mesh->nodes[i]->draw();

    mesh->draw();

    glLoadIdentity();
}


QVector3D GraphicWindow::getMousePos (int x, int y)
{
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;
    winX = float(x);
    winY = float(viewport[3]) - float(y);

    winZ = 1.;

    gluUnProject( winX, winY, winZ, model_view, projection, viewport,
                  &posX, &posY, &posZ);

    //        qWarning(" +++ GetMousePos  mouse position:  (%g, %g, %g) (%g, %g, %g)",
    //                 winX, winY, winZ, posX, posY, posZ);
    return QVector3D(posX, posY, posZ);
}


QVector2D GraphicWindow::normalizeMouse(QPoint qp)
{
    QVector2D res = QVector2D(arcb->invWidth * (double(qp.x()) - 0.5*screenW),
                              arcb->invHeight * (0.5*screenH - double(qp.y())));

    return res;
}



void GraphicWindow::euler2matr ()
{
    // Assume the angles are in passed in
    // in radians.
    double ch = cos (thetaRot*deg2rad); // heading
    double sh = sin (thetaRot*deg2rad);
    double ca = cos (psiRot  *deg2rad);   // attitude
    double sa = sin (psiRot  *deg2rad);
    double cb = cos (phiRot  *deg2rad);   // bank
    double sb = sin (phiRot  *deg2rad);

    transformM[0]  = ch*ca;
    transformM[1]  = sh*sb - ch*sa*cb;
    transformM[2]  = ch*sa*sb + sh*cb;
    transformM[3]  = 0.;
    transformM[4]  = sa;
    transformM[5]  = ca*cb;
    transformM[6]  = -ca*sb;
    transformM[7]  = 0.;
    transformM[8]  = -sh*ca;
    transformM[9]  = sh*sa*cb + ch*sb;
    transformM[10] = -sh*sa*sb + ch*cb;
    transformM[11] = 0.;
    transformM[12] = 0.;
    transformM[13] = 0.;
    transformM[14] = 0.;
    transformM[15] = 1.;
}

void GraphicWindow::quat2matr (QQuaternion q)
{
    double sqw = q.scalar()*q.scalar();
    double sqx = q.x()*q.x();
    double sqy = q.y()*q.y();
    double sqz = q.z()*q.z();

    double invs = 1 / (sqx + sqy + sqz + sqw);
    transformM[0] = (sqx - sqy - sqz + sqw)*invs;
    transformM[5] = (-sqx + sqy - sqz + sqw)*invs;
    transformM[10] = (-sqx - sqy + sqz + sqw)*invs;

    double tmp1 = q.x()*q.y();
    double tmp2 = q.z()*q.scalar();
    transformM[4] = 2.*(tmp1 + tmp2)*invs;
    transformM[1] = 2.*(tmp1 - tmp2)*invs;
    tmp1 = q.x()*q.z();
    tmp2 = q.y()*q.scalar();
    transformM[8] = 2.*(tmp1 - tmp2)*invs;
    transformM[2] = 2.*(tmp1 + tmp2)*invs;
    tmp1 = q.y()*q.z();
    tmp2 = q.x()*q.scalar();
    transformM[9] = 2.*(tmp1 + tmp2)*invs;
    transformM[6] = 2.*(tmp1 - tmp2)*invs;
    transformM[3] = transformM[7] = transformM[11] =
            transformM[12] = transformM[13] = transformM[14] = 0.;
    transformM[15] = 1.;
}


QVector3D GraphicWindow::quat2euler (QQuaternion q)
{
    double head, att, b;
    double test = q.x()*q.y() + q.z()*q.scalar();
    if (test > 0.499)
    { // singularity at north pole
        head = 2.*atan2(q.x(), q.scalar());
        att = pi_2;
        b = 0.;
        return QVector3D(head, att, b);
    }
    if (test < -0.499)
    { // singularity at south pole
        head = -2.*atan2(q.x(), q.scalar());
        att = - pi_2;
        b = 0.;
        return QVector3D(head, att, b);
    }
    double sqx = q.x()*q.x();
    double sqy = q.y()*q.y();
    double sqz = q.z()*q.z();
    head = atan2(2.*q.y()*q.scalar()-2.*q.x()*q.z(), 1. - 2.*sqy - 2.*sqz);
    att = asin(2.*test);
    b = atan2(2.*q.x()*q.scalar()-2.*q.y()*q.z(), 1. - 2.*sqx - 2.*sqz);
    return QVector3D(head, att, b);
}

QQuaternion GraphicWindow::quatfromEuler ()
{
    //   QQuaternion res;
    double heading2  = 0.5 * deg2rad * thetaRot;
    double attitude2 = 0.5 * deg2rad * psiRot;
    double bank2     = 0.5 * deg2rad * phiRot;

    double c1   = cos (heading2);
    double s1   = sin (heading2);
    double c2   = cos (attitude2);
    double s2   = sin (attitude2);
    double c3   = cos (bank2);
    double s3   = sin (bank2);
    double c1c2 = c1*c2;
    double s1s2 = s1*s2;
    return QQuaternion(c1c2*c3 - s1s2*s3,
                       c1c2*s3 + s1s2*c3, s1*c2*c3 + c1*s2*s3, c1*s2*c3 - s1*c2*s3);
}
