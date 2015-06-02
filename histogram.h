#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <QWidget>

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
class QHBoxLayout;
class QVBoxLayout;
QT_END_NAMESPACE


//! [0]
class Histogram : public QWidget
{
    Q_OBJECT

	int width;
	int height;
	bool setHistogram;

	int *histogrammArray;
	int *PhistogrammArray;

public:
    Histogram(QWidget *parent = 0, QImage *img = 0);
	void buildReadHistogramm();

public Q_SLOTS:
	void readHistogramm(QImage &img);

protected:
	void paintEvent(QPaintEvent *event);
private:


};
//! [0]



class Histograms : public QWidget
{
	Q_OBJECT


	QHBoxLayout *hbox;
	
public:
	Histogram hist1;
	Histogram hist2;
	
	Histograms(QWidget *parent = 0);
		

};


#endif
