#ifndef APPWINDOW_H
#define APPWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
class QHBoxLayout;
class QVBoxLayout;
QT_END_NAMESPACE

#include "imageviewer.h"
#include "histogram.h"

//! [0]
class AppWindow : public QMainWindow
{
    Q_OBJECT

public:
    AppWindow();

private:
	ImageViewer *view1;
	ImageViewer *view2;
	QHBoxLayout *hbox1;
	QHBoxLayout *hbox2;
	QVBoxLayout *vbox;

	Histogram *hist1;
	Histogram *hist2;
	//QWidget *hist;
	//QWidget *hist2;

};
//! [0]

#endif
