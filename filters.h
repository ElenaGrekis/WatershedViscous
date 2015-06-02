
#ifndef FILTERS_H
#define FILTERS_H

#include <QWidget>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace cv;

#include "graphs.h"

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
class QHBoxLayout;
class QVBoxLayout;
class QInputDialog;
QT_END_NAMESPACE


static IplImage* QImage2IplImage(QImage *qimg);

//! [0]
class Filters : public QWidget
{
	
    Q_OBJECT

	static const int m = 32;

	//Canny
	static  float GAUSSIAN_CUT_OFF;
	static  float MAGNITUDE_SCALE;
	static  float MAGNITUDE_LIMIT;
	static  int MAGNITUDE_MAX;
	static  int gaussianKernelWidth;
	//

	static int width;
	static int height;
	static int picsize;

	static int *data;
	 static int *magnitude;

	static float *xConv;
	static float *yConv;
	static float *xGradient;
	static float *yGradient;

public:
    Filters(QWidget *parent = 0);
	static double *svertka(double *val, int val_size, double *filt, int filt_size);
	static double *lp(double f0, const int M, double dt);
	static double *hp(double f0, const int M, double dt);
	static double* bp(double f1, double f2, const int M, double dt);
	static double* bs(double f1, double f2, const int M, double dt);

	static void main_filter(QImage &img, Filters2 fil, double fc, double fc2 = 0.0 );

	static void equilisation(QImage &img, bool ignoreBackground = false);

	static void dilate(QImage &img, int mask_size);
	static void erosion(QImage &img, int mask_size);
	static void DilDifEro(QImage &img, int dil_mask, int ero_mask);
	static void laplace(QImage &img);
	static void subOriginal(QImage &img1, QImage &orig);
	static int gradient (const int *data);
	static void gradientEdge(QImage &img);
	static void threadshold(QImage &img, int th);
	static void mopOpen(QImage &img, int dil_mask, int ero_mask);
	static void mopClose(QImage &img, int dil_mask, int ero_mask);

	static void minMax(const QImage &img, int &min, int &max);

	static void fillOne(int **mask, int w, int h);

	static int countObjects(QImage &img, int mask_size);
	
	static void gaussianBlur(QImage &img, float radius);
	

	static void visWatershed(QImage &img, QImage &markers, QString &);

	static void sobelThresh(QImage &pic);

	static void LoG(QImage &img, float blurRadius, float kernelRadius);

	static void zeroCrossing(QImage &img, float blurRadius, float kernelRadius);

	static void multiply(QImage &img1, QImage &orig);

	static void nullify(QImage &img1);

	static void combine(QImage &orig, QImage &img1);

	//Canny functions
	static void canny(QImage &img, float radius, float low, float high);
	static void computeGradients(float kernelRadius, int kernelWidth);
	static void performHysteresis(int low, int high);
	static void follow(int x1, int y1, int i1, int threshold);
	static void thresholdEdges();

 	static inline float gaussian(float x, float sigma)
	{
		return (float) exp(-(x * x) / (2.0f * sigma * sigma));
	}

	static inline float round(float number)
	{
		return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
	}
};
//! [0]

typedef struct CvWSNode
{
	struct CvWSNode* next;
	int mask_ofs;
	int img_ofs;
}
CvWSNode;

typedef struct CvWSQueue
{
	CvWSNode* first;
	CvWSNode* last;
}
CvWSQueue;


static CvWSNode* icvAllocWSNodes( CvMemStorage* storage );
#endif
