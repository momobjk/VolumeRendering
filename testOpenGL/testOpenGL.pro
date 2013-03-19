#-------------------------------------------------
#
# Project created by QtCreator 2013-03-04T11:42:23
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = testOpenGL
TEMPLATE = app


SOURCES += main.cpp \
    myWindow.cpp \
    myGLWidget.cpp

HEADERS  += \
    myWindow.h \
    myGLWidget.h

FORMS    += mainwindow.ui

LIBS = -lglut -lGLU
