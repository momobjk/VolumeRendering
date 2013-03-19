#include <QApplication>
#include "myWindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    myWindow myWin;
    myWin.show();
    return a.exec();
}
