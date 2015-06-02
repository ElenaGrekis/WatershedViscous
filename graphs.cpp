#include "graphs.h"
#include "canvas.h"

Graph::Graph(Filters2 filters, double *_file, int _filesize, QWidget *parent) : QWidget(parent)
{
	filter = filters;
	filter_dt = 0.001;

	arr = arr2 = result = nullptr;

	if(_file)
	{
		file = _file;
		size = _filesize;
	}
	else
	{
		file = nullptr;
		size = 0;
	}

	fft = new FFTRealWrapper; //fft for test

	QWidget *w = new QWidget( this );

	QGridLayout *hLayout = new QGridLayout;

	for(int i = 0; i < num; i++)
	{
		plot[i] = new Plot(this, graph_names[i]);
		//plot[i]->axisAutoScale(
		plot[i]->setAxisScale(QwtPlot::xBottom, 0, 1000);
		if(i == 1)
			plot[i]->setAxisScale(QwtPlot::yLeft, -0.9, 1.2);
		else
			plot[i]->setAxisScale(QwtPlot::yLeft, -100, 150);			

		QwtPlotGrid *grid = new QwtPlotGrid();
		grid->attach( plot[i] );

		curve[i] = new QwtPlotCurve("f(x)");

		curve[i]->setLegendAttribute( QwtPlotCurve::LegendShowLine, true );
		curve[i]->setStyle(QwtPlotCurve::Lines);
		curve[i]->setPen(QPen(Qt::blue));
		curve[i]->attach(plot[i]);    


		plot[i]->replot();

	}


	for(int i = 0, counter = 0; i < num; i++)
	{
		hLayout->addWidget( plot[i], counter, i & 1);
		if(i&1)
			counter+=2;
	}


	QHBoxLayout *horizontalButton = new QHBoxLayout;
	horizontalButton->addStretch(5);
	QPushButton *button = new QPushButton("Redraw");
	horizontalButton->addWidget(button);
	horizontalButton->addStretch(5);


	QVBoxLayout *mainVbox = new QVBoxLayout; //Main layout
	mainVbox->addLayout(hLayout);
	mainVbox->addLayout(horizontalButton);

	connect(button, SIGNAL(clicked()), this, SLOT(setData()));
	
	fBox1 = new QDoubleSpinBox;
	fBox1->setValue(20.0);
	fBox1->setPrefix("f1: ");
	fBox1->setMaximum(5000);
	horizontalButton->addWidget(fBox1);
	horizontalButton->addStretch(5);

	fBox2 = new QDoubleSpinBox;
	fBox2->setValue(70.0);
	fBox2->setPrefix("f2: ");
	fBox2->setMaximum(5000);
	horizontalButton->addWidget(fBox2);
	horizontalButton->addStretch(5);


	if(filter == HPFPOLY)
	{
		fBox1->setValue(80);
	}
	if(filter == BPFPOLY)
	{
		fBox1->setValue(10);
		fBox2->setValue(108);
	}
	if(filter == BSFPOLY)
	{
		fBox1->setValue(10);
		fBox2->setValue(108);
	}


	setData();

	setLayout(mainVbox);

//	setCentralWidget( w );
//	showMaximized();
}

double *svertka(double *val, int val_size, double *filt, int filt_size)
{
	
	int n = val_size,
		L = filt_size;
	double* arr_tmp = new double[n + L - 1];

	for(int i = 0; i < n + L - 1; i++)
		arr_tmp[i] = 0.0;

	for(int i = 0; i < n; i++)
		for(int j = 0; j < L; j++)
		{
			arr_tmp[i+j] += val[i] * filt[j];
		}
		/*
		int n = val_size,
		L = filt_size;
		double* arr_tmp = new double[n + L - 1];
		for(int i = 0; i < n + L - 1; i++)
		arr_tmp[i] = 0.0;
		int j0 = 0, jj = 0;

		for (int i = -L + 1; i < n; i++)
		{
			arr_tmp[i + L - 1] = 0;
			if (i > 0) j0 = i;
			for (int k = 0; k < L; k++)
			{
				jj = j0 + k;
				if (jj < i + L && jj < n)
				{
					arr_tmp[i + L - 1] += filt[L - 1 - k] * val[jj];
				}
				else break;
			}
		}

		double *res = new double [n];
		for (int i = 0; i < n; i++)
		res[i] = arr_tmp[i + 1 + L / 2]; 

		delete [] arr_tmp;
		
		return res;*/
		return &arr_tmp[L/2];
}

/*void MainWindow::conv_f(double re1, double im1, double re2, double im2, double &reOut, double &imOut)
{
reOut = re1 * re2 - im1 * im2;
imOut = -re1 * im2 - im1 * re2;
}*/
//#if 0
double * lp(double f0, const int M, double dt) 
{
	double d[] = {0.35577019, 0.2436983, 0.07211497, 0.00630165 };
	double *wlp = new double [2 * m + 1];// = new double;
	double tmp[m + 1];// = new double;
	double fact = 2 * f0 * dt;
	tmp[0] = fact;
	fact = fact * M_PI;
	for (int i = 1; i < m + 1; i++)
		tmp[i] = sin(fact * i) / M_PI / i;
	tmp[m] /= 2;

	double sumg = tmp[0], fact1, sum;
	for (int i = 1; i < m + 1; i++)
	{
		sum = d[0];
		fact1 = M_PI * i / m;
		for (int k = 1; k < 4; k++)
			sum += 2 * d[k] * cos(fact1 * k);
		tmp[i] *= sum;
		sumg += 2 * tmp[i];
	}
	for (int i = 1; i < m; i++)
		tmp[i] /= sumg;


	for (int i = 0; i < m; i++)
		wlp[m - 1 - i] = wlp[m + 1 + i] = tmp[i + 1];
	wlp[m] = tmp[0];
	return wlp;



}

double *hp(double f0, const int M, double dt)
{
	double *whp = new double [2 * m + 1];
	double *wlp = lp(f0, m, dt);
	for (int i = 1; i < m + 1; i++)
		whp[m + i] = whp[m - i] = - wlp[m + i];
	whp[m] = 1 - wlp[m];

	delete [] wlp;
	return whp;
}

double* bp(double f1, double f2, const int M, double dt)
{
	double* wbp = new double[2 * m + 1];
	double* wlp1 = lp(f1, m, dt);
	double* wlp2 = lp(f2, m, dt);
	for (int i = 0; i < m + 1; i++)
		wbp[m + i] = wbp[m - i] = wlp2[m + i] - wlp1[m + i];

	delete [] wlp1;
	delete [] wlp2;
	return wbp;
}
double* bs(double f1, double f2, const int M, double dt)
{
	double* wbs = new double[2 * m + 1];
	double* wlp1 = lp(f1, m, dt);
	double* wlp2 = lp(f2, m, dt);
	for (int i = 1; i < m + 1; i++)
		wbs[m + i] = wbs[m - i] = wlp1[m + i] - wlp2[m + i];
	wbs[m] = 1 + wlp1[m] - wlp2[m];

	delete [] wlp1;
	delete [] wlp2;
	return wbs;
}


double *Furie(double* arr_t, int arr_t_size, int max_f)
{
	int n = arr_t_size;
	double *res = new double[max_f];
	double re, im;
	double ddd = 1.0/22050;
	for (int ff = 0; ff < max_f; ff++)
	{
		re = im = 0;
		for (int i = 0; i < n; i++)
		{
			re += arr_t[i] * cos(M_PI * 2.0 * ff * i * ddd);
			im += arr_t[i] * sin(M_PI * 2.0 * ff * i * ddd);
		}
		res[ff] = sqrt(re * re + im * im) / n;
	}
	return res;
}

double *FillPolyharm(double X1, double X2, double X3, double f1, double f2, double f3, int n)//получение полигармонического процесса в массиве arrPolyharm[num]
{
	double *poly = new double[n];
	for (int i = 0; i < n; i++)
	{
		//if(i <= 160)
			poly[i] = X1 * sin(2.0 * M_PI * f1 * i * 0.001) +
			X2 * sin(2.0 * M_PI * f2 * i * 0.001) +
			X3 * sin(2.0 * M_PI * f3 * i * 0.001);
		//else poly[i] = 0.0;
	}
	return poly;
}

void Graph::setFile(double *file, int size)
{
	this->file = file;
	this->size = size;
	setData();
}

void Graph::setData()
{
	delete [] arr;
	delete [] arr2;
	delete [] result;

	arr = arr2 = result = nullptr;
	int count = 1000;

	double *arr4 = nullptr;
	double *arr3 = nullptr;


	double *data = nullptr;
	if(file)
	{
		data = file;
		count = size;
	}
	else
	{
		count = 1000;
		data = FillPolyharm(100.0, 20.0, 15.0, 75.0, 5.0, 160.0, count);
	}
	//double * data = FillPolyharm(10, 100, 15, 15, 60, 150, count);
	//double * data = FillPolyharm(15, 75, 20, 5, 50, 90, count);
	//int count = wavfile->read(m_buffer.data(), readLen);
	arr = new double[count];

	arr2 = data;

	for(int i = 0; i < count ; i++)
		arr[i] = i;

	if(filter != LPF && filter != HPF && filter != BPF && filter != BSF)
	{
		curve[0]->setSamples(arr, arr2, count );
		
		if(file)
		{
			plot[0]->setAxisAutoScale(0);
			plot[0]->setAxisAutoScale(1);
		}
		else
		{
		
			plot[0]->setAxisScale(QwtPlot::xBottom, 0, count);
			plot[0]->setAxisScale(QwtPlot::yLeft, -250, 400);
		}
		plot[0]->replot();
	}



	if(filter != LPF && filter != HPF && filter != BPF && filter != BSF)
	{
		double *result1 = new double[count];
		for(int i = 0; i < count; i++)
			result1[i] = 0.0;
		
		
		for(int i = 0; i < count; i++)
		{
			double res = 0.0, im = 0.0;
			fourier(arr2, i, res, im, count); 
			result1[i] = (sqrt(res*res + im*im) ) / count; //
		}

		//result1 = Furie(arr2, count, count);

		/*typedef FFTRealFixLenParam::DataType DataType;
			QVector<DataType> v1, v2;
			v1.resize(count);
			v2.resize(count);

			for(int i = 0; i < count; i++)
			{
				v1[i] = pcmToReal(arr2[i]);
			}

			fft->calculateFFT(v2.data(), v1.data()); //fft
			for(int i = 2; i <= count / 2; i++)
			{
				const qreal real = v2[i];
				qreal imag = 0.0;
				if (i > 0 && i < count / 2)
					imag = v2[count/2 + i];

				result1[i] = sqrt(real*real + imag*imag);//const qreal magnitude = sqrt(real*real + imag*imag);
			}
			*/

		curve[1]->setSamples(arr, result1, count / 2);
		if(file)
		{
			plot[1]->setAxisScale(QwtPlot::xBottom, 0, count / 2);
			plot[1]->setAxisScale(QwtPlot::yLeft, 0, 5);
			//plot[1]->setAxisAutoScale(QwtPlot::yLeft);
		}
		else
		{
			plot[1]->setAxisScale(QwtPlot::xBottom, 0, count / 2);
			plot[1]->setAxisScale(QwtPlot::yLeft, 0, 70);
		}
		plot[1]->replot();

		
		delete [] result1;
	}



	double f1 = fBox1->value();
	double f2 = fBox2->value();
	if(filter == LPF)
	{
		filter_dt = 0.005;
		arr3 = lp(f1, m, filter_dt);
	}
	else if(filter == HPF)
	{
		filter_dt = 0.005;
		arr3 = hp(f1, m, filter_dt);
	}
	else if(filter == BPF)
	{
		filter_dt = 0.005;
		arr3 = bp(f1, f2, m, filter_dt);
	}
	else if(filter == BSF)
	{
		filter_dt = 0.005;
		arr3 = bs(f1, f2, m, filter_dt);
	}
	else if(filter == LPFPOLY)
	{
		filter_dt = 0.005;
		arr3 = lp(f1, m, filter_dt);
	}
	else if(filter == HPFPOLY)
	{
		arr3 = hp(f1, m, filter_dt);
	}
	else if(filter == BPFPOLY)
		arr3 = bp(f1, f2, m, filter_dt);
	else if(filter == BSFPOLY)
		arr3 = bs(f1, f2, m, filter_dt);


	arr4 = svertka(arr2, count, arr3, 2*m + 1);


	if(filter == LPF || filter == HPF || filter == BPF || filter == BSF)
		curve[2]->setSamples(arr, arr3, 2*m); 
	else
		curve[2]->setSamples(arr, arr4, count); 

	if(filter == LPF || filter == HPF || filter == BPF || filter == BSF)
	{
		plot[2]->setAxisScale(QwtPlot::xBottom, 0, 2*m); 
		plot[2]->setAxisScale(QwtPlot::yLeft, -0.2, 1.0); 
	}
	else
	{
		if(file)
		{
			plot[2]->setAxisAutoScale(0);
			plot[2]->setAxisAutoScale(1);
		}
		else
		{
			plot[2]->setAxisScale(QwtPlot::xBottom, 0, count); 
			plot[2]->setAxisScale(QwtPlot::yLeft, -250, 400); 
		}
	}

	plot[2]->replot();

	delete [] result;
	result = nullptr;

	if(filter == LPF || filter == HPF || filter == BPF || filter == BSF)
	{
		result = new double[2*m];
		for(int i = 0; i < 2*m; i++)
		{
			double res = 0.0, im = 0.0;
			fourier(arr3, i, res, im, 2*m); 
			result[i] = sqrt(res*res + im*im) ;//
		}
		curve[3]->setSamples(arr, result, m);
		plot[3]->setAxisScale(QwtPlot::xBottom, 0, 2*m);
		plot[3]->setAxisScale(QwtPlot::yLeft, 0, 1.5);
	}
	else
	{
		result = new double[count];
		for(int i = 0; i < count; i++)
		{
			double res = 0.0, im = 0.0;
			fourier(arr4, i, res, im, count); 
			result[i] = sqrt(res*res + im*im) / count ;//
		}
		//result = Furie(arr4, count, count);

		curve[3]->setSamples(arr, result, count / 2);
		if(file)
		{
			plot[3]->setAxisScale(QwtPlot::xBottom, 0, count / 2);
			plot[3]->setAxisScale(QwtPlot::yLeft, 0, 5);
			//plot[3]->setAxisScale(QwtPlot::xBottom, 0, count / 2);
			//plot[3]->setAxisAutoScale(1);
		}
		else
		{
			plot[3]->setAxisScale(QwtPlot::xBottom, 0, count / 2);
			plot[3]->setAxisScale(QwtPlot::yLeft, 0, 70);
		}
	}

	//curve[3]->setRenderHint( QwtPlotItem::RenderAntialiased );

	plot[3]->replot();

	if(file)
		arr2 = nullptr;

	delete [] arr;
	delete [] arr2;
	delete [] arr3;
	//delete [] arr4;
	delete [] result;
	//delete [] data;

	arr = arr2 = arr3 = arr4 = result = nullptr;
}


double Graph::fourier(double *arr, double j, double &res, double &im, int arrsize)
{
	double r = 0, i = 0, am = 0, c = 0;


	for(int k = 0; k < arrsize; k++)
	{
		r += arr[k] * cos(2.0 * M_PI * k * j / arrsize);
		i += arr[k] * sin(2.0 * M_PI * k * j / arrsize);
	}
	res = r;
	im = i;


	return r - i;

}
