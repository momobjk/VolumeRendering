#include "myGLWidget.h"

myGLWidget::myGLWidget(int framesPerSecond, QWidget *parent, const std::string &name)
    : QGLWidget(parent)
{
    setWindowTitle(QString::fromStdString(name));
    if(framesPerSecond == 0)
        t_Timer = NULL;
    else
    {
        int seconde = 1000; // 1 seconde = 1000 ms
        int timerInterval = seconde / framesPerSecond;
        t_Timer = new QTimer(this);
        connect(t_Timer, SIGNAL(timeout()), this, SLOT(timeOutSlot()));
        t_Timer->start( timerInterval );        
    }

    b_Fullscreen = false;
}

void myGLWidget::keyPressEvent(QKeyEvent *keyEvent)
{
    switch(keyEvent->key())
    {
        case Qt::Key_Escape:
            close();
            break;
        case Qt::Key_F1:
            toggleFullWindow();
            break;
    }
}

void myGLWidget::timeOutSlot()
{
}

void myGLWidget::toggleFullWindow()
{
    if(b_Fullscreen)
    {
        showNormal();
        b_Fullscreen = false;
    }
    else
    {
        showFullScreen();
        b_Fullscreen = true;
    }
}
