#include "mainwindow.h"
//#include "appwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
	tab = new QTabWidget(this);

	tab->addTab( new Graph(LPF), "LPF" );
	tab->addTab( new Graph(HPF), "HPF" );
	tab->addTab( new Graph(BPF), "BPF" );
	tab->addTab( new Graph(BSF), "BSF" );
	
	Images *img = new Images;
	Histograms *hist = new Histograms;
	Graph *gr = new Graph(LPFPOLY);


	tab->addTab( img, "Image" );
	tab->addTab( hist, "Histograms" );


	tab->setCurrentIndex(4);
	
	Q_ASSERT(connect(img->view1, SIGNAL(sendImage(QImage &)),
		&hist->hist1, SLOT(readHistogramm(QImage &)) ));

	Q_ASSERT(connect(img->view2, SIGNAL(sendImage(QImage &)),
		&hist->hist2, SLOT(readHistogramm(QImage &)) ));

	Q_ASSERT(connect(img->view1, SIGNAL(sendImage(QImage &)),
		img->view2, SLOT(setImage(QImage &)) ));

	Q_ASSERT(connect(img->view2, SIGNAL(updateImageA(QImage &)),
		img->view1, SLOT(setImage(QImage &)) ));

	Q_ASSERT(connect(img->buttonSNPNoise, SIGNAL(clicked()),
		img->view2, SLOT(addSNPNoise()) ));

	Q_ASSERT(connect(img->buttonGaussNoise, SIGNAL(clicked()),
		img->view2, SLOT(addGaussNoise()) ));

	Q_ASSERT(connect(img->buttonLPF, SIGNAL(clicked()),
		img->view2, SLOT(lpfFilter()) ));

	Q_ASSERT(connect(img->buttonHPF, SIGNAL(clicked()),
		img->view2, SLOT(hpfFilter()) ));

	Q_ASSERT(connect(img->buttonBPF, SIGNAL(clicked()),
		img->view2, SLOT(bpfFilter()) ));

	Q_ASSERT(connect(img->buttonBSF, SIGNAL(clicked()),
		img->view2, SLOT(bsfFilter()) ));

	Q_ASSERT(connect(img->equlise, SIGNAL(clicked()),
		img->view2, SLOT(equilise()) ));

	Q_ASSERT(connect(img->dilate, SIGNAL(clicked()),
		img->view2, SLOT(dilate()) ));

	Q_ASSERT(connect(img->erosion, SIGNAL(clicked()),
		img->view2, SLOT(erosion()) ));

	Q_ASSERT(connect(img->dilDifEro, SIGNAL(clicked()),
		img->view2, SLOT(DilDifEro()) ));

	Q_ASSERT(connect(img->laplace, SIGNAL(clicked()),
		img->view2, SLOT(laplace()) ));

	Q_ASSERT(connect(img->view2, SIGNAL(getOriginalImage(QImage &)),
		img->view1, SLOT(putOriginalImage(QImage &)) ));

	Q_ASSERT(connect(img->subBfromA, SIGNAL(clicked()),
		img->view2, SLOT(subBfromA()) ));
	
	Q_ASSERT(connect(img->subAfromB, SIGNAL(clicked()),
		img->view2, SLOT(subAfromB()) ));

	Q_ASSERT(connect(img->gradientFilter, SIGNAL(clicked()),
		img->view2, SLOT(gradientFilter()) ));

	Q_ASSERT(connect(img->threadshold, SIGNAL(clicked()),
		img->view2, SLOT(threadshold()) ));

	Q_ASSERT(connect(img->morphOpen, SIGNAL(clicked()),
		img->view2, SLOT(mopOpen()) ));

	Q_ASSERT(connect(img->morphClose, SIGNAL(clicked()),
		img->view2, SLOT(mopClose()) ));

	Q_ASSERT(connect(img->copyAB, SIGNAL(clicked()),
		img->view2, SLOT(copyAB()) ));

	Q_ASSERT(connect(img->copyBA, SIGNAL(clicked()),
		img->view2, SLOT(copyBA()) ));


	Q_ASSERT(connect(img->saveB, SIGNAL(clicked()),
		img->view2, SLOT(saveB()) ));

	Q_ASSERT(connect(img->count, SIGNAL(clicked()),
		img->view2, SLOT(countObjects()) ));

	Q_ASSERT(connect(img->gaussianBlur, SIGNAL(clicked()),
		img->view2, SLOT(gaussianBlur()) ));

	Q_ASSERT(connect(img->visWater, SIGNAL(clicked()),
		img->view2, SLOT(visWatershed()) ));

	Q_ASSERT(connect(img->sobelBtn, SIGNAL(clicked()),
		img->view2, SLOT(sobel()) ));

	Q_ASSERT(connect(img->cannyBtn, SIGNAL(clicked()),
		img->view2, SLOT(canny()) ));

	Q_ASSERT(connect(img->invertBtn, SIGNAL(clicked()),
		img->view2, SLOT(invert()) ));

	Q_ASSERT(connect(img->logBtn, SIGNAL(clicked()),
		img->view2, SLOT(LoG()) ));

	Q_ASSERT(connect(img->zeroCrossBtn, SIGNAL(clicked()),
		img->view2, SLOT(zeroCrossing()) ));
	
	Q_ASSERT(connect(img->multBtn, SIGNAL(clicked()),
		img->view2, SLOT(multiply()) ));

	Q_ASSERT(connect(img->nullifyBtn, SIGNAL(clicked()),
		img->view2, SLOT(nullify()) ));
	
	Q_ASSERT(connect(img->combineBtn, SIGNAL(clicked()),
		img->view2, SLOT(combine()) ));

	setCentralWidget( tab );
	showMaximized();
}

