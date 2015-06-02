######################################################################
# Automatically generated by qmake (2.01a) ?? 28. ??? 00:28:28 2014
######################################################################

CONFIG += qwt

TEMPLATE = app
TARGET = imageviewer4
DEPENDPATH += .
INCLUDEPATH += .

QT += multimedia

CONFIG += ordered

SUBDIRS += fftreal

INCLUDEPATH += fftreal
LIBS += -L.. -lfftreal

# Input
HEADERS += filters.h noise.h canvas.h graphs.h histogram.h imageviewer.h mainwindow.h
SOURCES += filters.cpp noise.cpp \
           canvas.cpp \
           graphs.cpp \
           histogram.cpp \
           imageviewer.cpp \
           main.cpp \
           mainwindow.cpp