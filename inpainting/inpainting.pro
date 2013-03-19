#-------------------------------------------------
#
# Project created by QtCreator 2013-03-07T14:48:14
#
#-------------------------------------------------

QMAKE_CXXFLAGS += -std=c++11

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = inpainting
TEMPLATE = app


SOURCES += main.cpp \
    view.cpp \
    glWidget.cpp

HEADERS  += \
    view.h \
    glWidget.h

FORMS    += \
    view.ui

LIBS = -lglut -lGLU

