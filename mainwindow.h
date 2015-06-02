#pragma once

#include <QMainWindow>
#include <QWidget>
#include <QtGui>
#include <QTabWidget>
#include "graphs.h"
#include "histogram.h"
#include "imageviewer.h"
#include "noise.h"


class MainWindow : public QMainWindow
{
	Q_OBJECT

		
	QTabWidget *tab;

public:
	MainWindow(QWidget *parent = 0);
	double *readFile(int &);
	double *GetValuesFromFile(char *arr, int &size);





};