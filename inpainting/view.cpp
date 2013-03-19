#include "view.h"


View::View(QWidget* parent, Qt::WindowFlags flags): QMainWindow(parent, flags)
{
  setupUi(this);
  setupActions();
  glWindow = new GLWidget(60, widget);
}

View::~View()
{

}

void View::setupActions()
{
  connect(actionOpen,SIGNAL(triggered()), this, SLOT(open()));
  connect(pushButtonOK ,SIGNAL(clicked()), this, SLOT(printSizes()));
}


void View::open()
{  
   QString fileName = QFileDialog::getOpenFileName(this,
                                    tr("Open File"), QDir::currentPath());
   std::cout << fileName.toStdString() << std::endl;
    if (!fileName.isEmpty()) {
        QImage image(fileName);
        if (image.isNull()) {
            QMessageBox::information(this, tr("Image Viewer"),
                                     tr("Cannot load %1.").arg(fileName));
            return;
    }
	      
// 	  scrollArea->adjustSize();
  }
}


void View::keyPressEvent(QKeyEvent *keyEvent)
{
    switch(keyEvent->key())
    {
        case Qt::Key_Escape:
            close();
            break;
    }
}

void View::resizeEvent(QResizeEvent *resEvent){
    if(resEvent != NULL)
    {
        widget->setMinimumHeight(centralwidget->height() - 30);
        widget->setMaximumHeight(centralwidget->height() - 30);

        glWindow->setMinimumHeight(widget->height());
        glWindow->setMinimumWidth(widget->width());
        glWindow->setMaximumHeight(widget->height());
        glWindow->setMaximumWidth(widget->width());
    }
}

void View::printSizes()
{
    std::cout<<"MainWindow : " +  std::to_string(this->width())  + " x " +  std::to_string(this->height()) << std::endl;
    std::cout<<"centralWidget : " +  std::to_string(centralwidget->width())  + " x " +  std::to_string(centralwidget->height()) << std::endl;
    std::cout<<"widget : " +  std::to_string(widget->width())  + " x " +  std::to_string(widget->height()) << std::endl;
    std::cout<<"glWindow : " + std::to_string(glWindow->width())  + " x " + std::to_string(glWindow->height()) << std::endl;

}
