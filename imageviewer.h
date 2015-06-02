#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QMainWindow>
#include <ctime>

#include "noise.h"
#include "filters.h"

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
class QHBoxLayout;
class QVBoxLayout;
class QPushButton;
class QSpinBox;
class QDoubleSpinBox;
class QGroupBox;
QT_END_NAMESPACE

class Histogram;
class Images;
//! [0]
class ImageViewer : public QMainWindow
{
    Q_OBJECT

	Images *parent;

Q_SIGNALS:
	void sendImage(QImage &img);
	void fileOpened(double *file, int size);
	void getOriginalImage(QImage &img);
	void updateImageA(QImage &img);
public:
    ImageViewer(QWidget * parent = 0, Histogram * = 0);
	double* convertImageToDouble(QImage &img);
	QImage convertDoubleToImage(double *img, int width, int height);

public Q_SLOTS:
	void putOriginalImage(QImage &img);
	void setImage(QImage &img);
	void addSNPNoise();
	void addGaussNoise();
	void lpfFilter();
	void hpfFilter();
	void bpfFilter();
	void bsfFilter();
	void equilise();
	void dilate();
	void erosion();
	void DilDifEro();
	void laplace();
	void subBfromA();
	void subAfromB();
	void copyBA();
	void copyAB();
	void saveB();
	void gradientFilter();
	void threadshold();
	void mopOpen();
	void mopClose();
	void countObjects();
	void gaussianBlur();
	void visWatershed();
	void sobel();
	void canny();
	void invert();
	void LoG();
	void zeroCrossing();
	void multiply();
	void nullify();
	void combine();
	
private Q_SLOTS:
    void open();
    void zoomIn();
    void zoomOut();
    void normalSize();
    void fitToWindow();

private:
    void createActions();
    void createMenus();
    void updateActions();
    void scaleImage(double factor);
    void adjustScrollBar(QScrollBar *scrollBar, double factor);

    QLabel *imageLabel;
    QScrollArea *scrollArea;
    double scaleFactor;

    QAction *openAct;
    QAction *zoomInAct;
    QAction *zoomOutAct;
    QAction *normalSizeAct;
    QAction *fitToWindowAct;

    QMenu *fileMenu;
    QMenu *viewMenu;
	
};
//! [0]

class Images : public QWidget
{

	QHBoxLayout *hbox1;

public:
	ImageViewer *view1;
	ImageViewer *view2;

	QPushButton *buttonSNPNoise;
	QPushButton *buttonGaussNoise;
	QDoubleSpinBox *snpPerc;
	QDoubleSpinBox *gaussRazb;
	QDoubleSpinBox *gaussianRadius;
	QSpinBox *dilateMask;
	QSpinBox *eroMask;
	QSpinBox *objMask;

	QPushButton *buttonLPF;
	QPushButton *buttonHPF;
	QPushButton *buttonBPF;
	QPushButton *buttonBSF;
	QPushButton *equlise;
	QCheckBox *equliseBox;
	QPushButton *dilate;
	QPushButton *erosion;
	QPushButton *dilDifEro;
	QPushButton *laplace;
	QPushButton *subBfromA;
	QPushButton *subAfromB;
	QPushButton *gradientFilter;
	QDoubleSpinBox *spinboxFreq;
	QDoubleSpinBox *spinboxFreq2;
	QPushButton *threadshold;
	QSpinBox *threadsholdValue;
	QPushButton *morphOpen;
	QPushButton *morphClose;
	QPushButton *copyAB;
	QPushButton *copyBA;
	QPushButton *saveB;
	QPushButton *count;
	QPushButton *gaussianBlur;
	QPushButton *visWater;
	QPushButton *sobelBtn;
	QPushButton *invertBtn;
	QPushButton *logBtn;
	QDoubleSpinBox *logBlurRadius;
	QSpinBox *logKernelRadius;
	QPushButton *zeroCrossBtn;
	QPushButton *multBtn;
	QPushButton *nullifyBtn;
	QPushButton *combineBtn;

	QPushButton *cannyBtn;
	QDoubleSpinBox *cannyRadius;
	QDoubleSpinBox *cannyLow;
	QDoubleSpinBox *cannyHigh;

	Images(QWidget * parent = 0);
};

Q_GUI_EXPORT void qt_blurImage( QImage &, qreal radius, bool quality, int transposed = 0);
#endif
