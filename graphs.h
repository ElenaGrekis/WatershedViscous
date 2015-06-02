#pragma once

#include <QMainWindow>
#include <QWidget>
#include <QtGui>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_curve_fitter.h>
#include <QGlobal.h>
#include <random>
#include <ctime>
#include "fftreal_wrapper.h"
#include "FFTRealFixLenParam.h"


const quint16 PCMS16MaxAmplitude =  32768; // because minimum is -32768

static const char* graph_names[] = 
{
	"",
	"Преобразование Фурье",
	"",
	""
};
static const int m = 32;
//static const int N = 1000;
//static const int M = 200;
static const int num = 4;
static const double global_dt = 0.001;//0.005;

enum Filters2 {LPF, HPF, BPF, BSF, LPFPOLY, HPFPOLY, BPFPOLY, BSFPOLY, WAV};

class Graph : public QWidget
{
	Q_OBJECT

	double alpha;
	double f;
	double strength;

	QwtPlot *plot[num];
    QwtPlotCurve *curve[num];

	double *arr;
	double *arr2;
	double *result;

	double *file;
	int size;

	Filters2 filter;
	double filter_dt;//0.005;
	QDoubleSpinBox *fBox1;
	QDoubleSpinBox *fBox2;

	FFTRealWrapper *fft;


	double fourier(double *arr, double j, double &res, double &im, int);

public:
	Graph(enum Filters2, double *file = nullptr, int size = 0, QWidget *parent = 0);

	public Q_SLOTS:
		void setData();
		void setFile(double *file, int size);

};