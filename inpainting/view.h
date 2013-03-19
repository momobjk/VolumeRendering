#ifndef VIEW_H
#define VIEW_H

#include "ui_view.h"


#include <QMainWindow>
#include<QWidget>
#include <QEvent>
#include <QMouseEvent>
#include <QLabel>
#include <QString>
#include <iostream>
#include <list>
#include <memory>
#include <QFileDialog>
#include <QDir>
#include <QMessageBox>

#include "glWidget.h"



class View : public QMainWindow, private Ui_MainWindow {
  Q_OBJECT
  
public:    
    View(QWidget* parent = 0, Qt::WindowFlags flags = 0);
    void keyPressEvent( QKeyEvent *keyEvent);
    void resizeEvent(QResizeEvent *resEvent);
    void setupActions();
    virtual ~View();

private:
  GLWidget* glWindow;

private slots:
  void open();
  void printSizes();
};

#endif
