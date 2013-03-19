#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QtOpenGL>
#include <QGLWidget>

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(int framesPerSecond = 0, QWidget *parent = 0);
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void toggleFullWindow();

public slots:
    virtual void timeOutSlot();

private:
    QTimer *t_Timer;
};


#endif // MYGLWIDGET_H
