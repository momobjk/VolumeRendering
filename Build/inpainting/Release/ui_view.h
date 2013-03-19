/********************************************************************************
** Form generated from reading UI file 'view.ui'
**
** Created: Tue Mar 19 13:41:18 2013
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_VIEW_H
#define UI_VIEW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionNew_Selection;
    QWidget *centralwidget;
    QGridLayout *gridLayout_5;
    QSpacerItem *verticalSpacer_3;
    QPushButton *pushButtonOK;
    QSpacerItem *horizontalSpacer_8;
    QLabel *labelWidth;
    QPushButton *pushButtonColor;
    QSpacerItem *horizontalSpacer_4;
    QSpacerItem *verticalSpacer_2;
    QSpacerItem *horizontalSpacer_5;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QSpacerItem *horizontalSpacer_2;
    QSpacerItem *horizontalSpacer_9;
    QSpinBox *spinBoxWidth;
    QSpacerItem *horizontalSpacer_3;
    QSpacerItem *horizontalSpacer_6;
    QSpacerItem *horizontalSpacer_7;
    QWidget *widget;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuSelection;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(936, 665);
        MainWindow->setMinimumSize(QSize(936, 665));
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionNew_Selection = new QAction(MainWindow);
        actionNew_Selection->setObjectName(QString::fromUtf8("actionNew_Selection"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        centralwidget->setMinimumSize(QSize(936, 620));
        gridLayout_5 = new QGridLayout(centralwidget);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        verticalSpacer_3 = new QSpacerItem(82, 425, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_5->addItem(verticalSpacer_3, 5, 3, 3, 1);

        pushButtonOK = new QPushButton(centralwidget);
        pushButtonOK->setObjectName(QString::fromUtf8("pushButtonOK"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(pushButtonOK->sizePolicy().hasHeightForWidth());
        pushButtonOK->setSizePolicy(sizePolicy);

        gridLayout_5->addWidget(pushButtonOK, 4, 3, 1, 1);

        horizontalSpacer_8 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_8, 3, 4, 1, 1);

        labelWidth = new QLabel(centralwidget);
        labelWidth->setObjectName(QString::fromUtf8("labelWidth"));
        sizePolicy.setHeightForWidth(labelWidth->sizePolicy().hasHeightForWidth());
        labelWidth->setSizePolicy(sizePolicy);
        QFont font;
        font.setPointSize(11);
        font.setBold(true);
        font.setWeight(75);
        labelWidth->setFont(font);

        gridLayout_5->addWidget(labelWidth, 1, 3, 1, 1);

        pushButtonColor = new QPushButton(centralwidget);
        pushButtonColor->setObjectName(QString::fromUtf8("pushButtonColor"));
        sizePolicy.setHeightForWidth(pushButtonColor->sizePolicy().hasHeightForWidth());
        pushButtonColor->setSizePolicy(sizePolicy);

        gridLayout_5->addWidget(pushButtonColor, 3, 3, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_4, 3, 2, 1, 1);

        verticalSpacer_2 = new QSpacerItem(5, 3, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_5->addItem(verticalSpacer_2, 0, 1, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_5, 4, 2, 1, 1);

        horizontalSpacer = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer, 5, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        gridLayout_5->addItem(verticalSpacer, 7, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_2, 1, 2, 1, 1);

        horizontalSpacer_9 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_9, 4, 4, 1, 1);

        spinBoxWidth = new QSpinBox(centralwidget);
        spinBoxWidth->setObjectName(QString::fromUtf8("spinBoxWidth"));
        sizePolicy.setHeightForWidth(spinBoxWidth->sizePolicy().hasHeightForWidth());
        spinBoxWidth->setSizePolicy(sizePolicy);
        spinBoxWidth->setMinimum(3);
        spinBoxWidth->setMaximum(20);

        gridLayout_5->addWidget(spinBoxWidth, 2, 3, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_3, 2, 2, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_6, 1, 4, 1, 1);

        horizontalSpacer_7 = new QSpacerItem(5, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_7, 2, 4, 1, 1);

        widget = new QWidget(centralwidget);
        widget->setObjectName(QString::fromUtf8("widget"));
        QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(widget->sizePolicy().hasHeightForWidth());
        widget->setSizePolicy(sizePolicy1);
        widget->setMinimumSize(QSize(812, 570));
        widget->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_5->addWidget(widget, 1, 1, 6, 1);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 936, 29));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuSelection = new QMenu(menubar);
        menuSelection->setObjectName(QString::fromUtf8("menuSelection"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuSelection->menuAction());
        menuFile->addAction(actionOpen);
        menuSelection->addAction(actionNew_Selection);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionNew_Selection->setText(QApplication::translate("MainWindow", "New Selection", 0, QApplication::UnicodeUTF8));
        pushButtonOK->setText(QApplication::translate("MainWindow", "ok", 0, QApplication::UnicodeUTF8));
        labelWidth->setText(QApplication::translate("MainWindow", "Width", 0, QApplication::UnicodeUTF8));
        pushButtonColor->setText(QApplication::translate("MainWindow", "Color", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuSelection->setTitle(QApplication::translate("MainWindow", "Selection", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_VIEW_H
