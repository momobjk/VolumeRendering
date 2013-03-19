#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QtOpenGL>
#include <QGLWidget>

class myGLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit myGLWidget(int framesPerSecond = 0, QWidget *parent = 0, const std::string &name = "");
    virtual void initializeGL() = 0;
    virtual void resizeGL(int width, int height) = 0;
    virtual void paintGL() = 0;
    virtual void keyPressEvent( QKeyEvent *keyEvent );
    void toggleFullWindow();

public slots:
    virtual void timeOutSlot();

private:
    QTimer *t_Timer;
    bool b_Fullscreen;

};


#endif // MYGLWIDGET_H
