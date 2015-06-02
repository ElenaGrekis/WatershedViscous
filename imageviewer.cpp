#include <QtGui>

#include "imageviewer.h"
#include "histogram.h"
#include "noise.h"


//! [0]
ImageViewer::ImageViewer(QWidget *parent, Histogram *_hist)
{
	this->parent = static_cast<Images*>(parent);

	imageLabel = new QLabel;
	imageLabel->setBackgroundRole(QPalette::Base);
	imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageLabel->setScaledContents(true);

	scrollArea = new QScrollArea;
	scrollArea->setBackgroundRole(QPalette::Dark);
	scrollArea->setWidget(imageLabel);
	setCentralWidget(scrollArea);

	createActions();
	createMenus();



	setWindowTitle(tr("Image Viewer"));
	resize(500, 400);
}
//! [0]

//! [1]
void ImageViewer::open()
	//! [1] //! [2]
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open File"), QDir::currentPath());
	if (!fileName.isEmpty()) 
	{

		QFileInfo info(fileName);
		QString x = info.suffix().toLower();
		if(info.suffix().toLower() == "img")
		{
			QByteArray imgData;
			QFile file(fileName);

			if (!file.open(QIODevice::ReadOnly))
			{
				QMessageBox::information(this, tr("Image Viewer"),
					tr("Cannot load %1.").arg(fileName));
				return;
			}

			int imgSizeHeight = QInputDialog::getInteger( this, "Height" , "Enter height: ", 1024, 0 );
			int imgSizeWidth = QInputDialog::getInteger( this, "Width" , "Enter width: ", 1024, 0 );

			if(!imgSizeHeight || !imgSizeWidth)
				return;

			file.seek(2048);
			imgData = file.readAll();

			unsigned long long q = 0;
			unsigned short min = USHRT_MAX;
			unsigned short max = 0;

			char* data = imgData.data();
			for(int i = 0 ; i < imgSizeHeight; i++)
				for(int j = 0; j < imgSizeWidth; j++)
				{

					unsigned short w = (((data[q] & 0x00ff) << 8)  )  | ((data[q+1] & 0x00ff ) );


					if(w < min)
						min = w;
					if(w > max)
						max = w;

					q += 2;
				}

				int x = 0;
				unsigned char *tmp = new unsigned char[imgSizeHeight*imgSizeWidth];

				unsigned char maxChar = 0;
				unsigned char minChar = 255;
				q = 0;
				for(int i = 0; i < imgSizeHeight; i++)
					for(int j = 0; j < imgSizeWidth; j++)
					{
						unsigned short w = (((data[q] & 0x00ff) << 8)  )  | ((data[q+1] & 0x00ff ) );

						double m = 255.0/ (max-min);
						double fl = (w * m) - (m*min) ;
						int n = fl;//((double)w / 0xffff) * 255.0; //fl



						if(n > 255)
							n = 255;
						else if(n < 0)
							n = 0; 

						if(n > maxChar)
							maxChar = n;
						if(n < minChar)
							minChar = n;
						/*if(i == 100 && j == 100)
							int pixel = tmp[i*imgSizeWidth + j];*/
						tmp[i*imgSizeWidth + j] = (static_cast<unsigned char>(n)) & 0x00ff;

						q += 2;
					}

					//инверсия картинки
					for(int i = 0; i < imgSizeHeight; i++)
						for(int j = 0; j < imgSizeWidth; j++)
						{
							tmp[i*imgSizeWidth + j] = maxChar - tmp[i*imgSizeWidth + j];
						}

						QImage image((unsigned char*)tmp, imgSizeHeight, imgSizeWidth, QImage::Format_Indexed8);

						QVector<QRgb> palette;
						for(int i=0; i<256; i++)
							palette.append(qRgb(i,i,i));

						image.setColorTable(palette);
						setImage(image);


						delete [] tmp;

		}
		else
		{
			QImage image(fileName);
			if (image.isNull()) {
				QMessageBox::information(this, tr("Image Viewer"),
					tr("Cannot load %1.").arg(fileName));
				return;
			}
			/*
			QImage::Format f = image.format();
			int maxChar = 0;
			int imgSizeHeight = image.height();
			int imgSizeWidth=image.width();
			for(int i = 0;i < imgSizeHeight; i++)
				for(int j = 0;j < imgSizeWidth; j++)
				{

					unsigned short q = QColor(image.pixel(i, j)).red();
					if(q > maxChar)
					{
						maxChar = q;
					}
				}
				*/
				//	int pixel = QColor(image.pixel(100, 100)).red();

				setImage(image);
		}

	}
}

double* ImageViewer::convertImageToDouble(QImage &img)
{
	int width = img.width();
	int height = img.height();
	double *result_img = new double[width*height];
	int count = 0;
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			/*img->pix[j][i].r = QColor(temp->pixel(i, j)).red();
			img->pix[j][i].g = QColor(temp->pixel(i, j)).red(); //green
			img->pix[j][i].b = QColor(temp->pixel(i, j)).red(); //blue*/
			result_img[count] = QColor(img.pixel(i, j)).red();
			count++;
		}

		return result_img;
}

QImage ImageViewer::convertDoubleToImage(double *img_result, int width, int height)
{
	QImage img(imageLabel->pixmap()->toImage()); //FIXME
	int count = 0;
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			img.setPixel(i, j, QColor(img_result[count], img_result[count], img_result[count]).rgba());
			count++;
		}
		return img;
}

void ImageViewer::setImage(QImage &img)
{
	imageLabel->setPixmap(QPixmap::fromImage(img));
	emit sendImage(img);

	//double *res = convertImageToDouble(img);
	//emit fileOpened(res, img.width() * img.height() );

	//! [3] //! [4]
	scaleFactor = 1.0;

	fitToWindowAct->setEnabled(true);
	updateActions();

	if (!fitToWindowAct->isChecked())
		imageLabel->adjustSize();
}

//! [9]
void ImageViewer::zoomIn()
	//! [9] //! [10]
{
	scaleImage(1.25);
}

void ImageViewer::zoomOut()
{
	scaleImage(0.8);
}

//! [10] //! [11]
void ImageViewer::normalSize()
	//! [11] //! [12]
{
	imageLabel->adjustSize();
	scaleFactor = 1.0;
}
//! [12]

//! [13]
void ImageViewer::fitToWindow()
	//! [13] //! [14]
{
	bool fitToWindow = fitToWindowAct->isChecked();
	scrollArea->setWidgetResizable(fitToWindow);
	if (!fitToWindow) {
		normalSize();
	}
	updateActions();
}
//! [14]

//! [17]
void ImageViewer::createActions()
	//! [17] //! [18]
{
	openAct = new QAction(tr("&Open..."), this);
	openAct->setShortcut(tr("Ctrl+O"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	zoomInAct = new QAction(tr("Zoom &In (25%)"), this);
	zoomInAct->setShortcut(tr("Ctrl++"));
	zoomInAct->setEnabled(false);
	connect(zoomInAct, SIGNAL(triggered()), this, SLOT(zoomIn()));

	zoomOutAct = new QAction(tr("Zoom &Out (25%)"), this);
	zoomOutAct->setShortcut(tr("Ctrl+-"));
	zoomOutAct->setEnabled(false);
	connect(zoomOutAct, SIGNAL(triggered()), this, SLOT(zoomOut()));

	normalSizeAct = new QAction(tr("&Normal Size"), this);
	normalSizeAct->setShortcut(tr("Ctrl+S"));
	normalSizeAct->setEnabled(false);
	connect(normalSizeAct, SIGNAL(triggered()), this, SLOT(normalSize()));

	fitToWindowAct = new QAction(tr("&Fit to Window"), this);
	fitToWindowAct->setEnabled(false);
	fitToWindowAct->setCheckable(true);
	fitToWindowAct->setShortcut(tr("Ctrl+F"));
	connect(fitToWindowAct, SIGNAL(triggered()), this, SLOT(fitToWindow()));

}
//! [18]

//! [19]
void ImageViewer::createMenus()
	//! [19] //! [20]
{
	fileMenu = new QMenu(tr("&File"), this);
	fileMenu->addAction(openAct);

	viewMenu = new QMenu(tr("&View"), this);
	viewMenu->addAction(zoomInAct);
	viewMenu->addAction(zoomOutAct);
	viewMenu->addAction(normalSizeAct);
	viewMenu->addSeparator();
	viewMenu->addAction(fitToWindowAct);

	menuBar()->addMenu(fileMenu);
	menuBar()->addMenu(viewMenu);
}
//! [20]

//! [21]
void ImageViewer::updateActions()
	//! [21] //! [22]
{
	zoomInAct->setEnabled(!fitToWindowAct->isChecked());
	zoomOutAct->setEnabled(!fitToWindowAct->isChecked());
	normalSizeAct->setEnabled(!fitToWindowAct->isChecked());
}
//! [22]

//! [23]
void ImageViewer::scaleImage(double factor)
	//! [23] //! [24]
{
	Q_ASSERT(imageLabel->pixmap());
	scaleFactor *= factor;
	imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());

	adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
	adjustScrollBar(scrollArea->verticalScrollBar(), factor);

	zoomInAct->setEnabled(scaleFactor < 3.0);
	zoomOutAct->setEnabled(scaleFactor > 0.333);
}
//! [24]

//! [25]
void ImageViewer::adjustScrollBar(QScrollBar *scrollBar, double factor)
	//! [25] //! [26]
{
	scrollBar->setValue(int(factor * scrollBar->value()
		+ ((factor - 1) * scrollBar->pageStep()/2)));
}
//! [26]

void ImageViewer::addSNPNoise()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Noise::saltPepper(tmp, parent->snpPerc->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::addGaussNoise()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Noise::gaussNoise(tmp, parent->gaussRazb->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::lpfFilter()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	Filters::main_filter(tmp, LPF, parent->spinboxFreq->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::hpfFilter()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	Filters::main_filter(tmp, HPF, parent->spinboxFreq->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::bpfFilter()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	Filters::main_filter(tmp, BPF, parent->spinboxFreq->value(), parent->spinboxFreq2->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::bsfFilter()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	Filters::main_filter(tmp, BSF, parent->spinboxFreq->value(), parent->spinboxFreq2->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::equilise()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::equilisation(tmp, parent->equliseBox->isChecked());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::dilate()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::dilate(tmp, parent->dilateMask->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::erosion()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::erosion(tmp, parent->eroMask->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::DilDifEro()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::DilDifEro(tmp, parent->dilateMask->value(), parent->eroMask->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::laplace()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::laplace(tmp);

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::putOriginalImage(QImage &img)
{
	img = imageLabel->pixmap()->toImage();
}

void ImageViewer::subBfromA()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage original(tmp);

	emit getOriginalImage(original);


	Filters::subOriginal(tmp, original);

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}


void ImageViewer::subAfromB()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage original(tmp);

	emit getOriginalImage(original);


	Filters::subOriginal(original, tmp);

	imageLabel->setPixmap(QPixmap::fromImage(original)); //FIXME

	sendImage(tmp);
}

void ImageViewer::threadshold()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::threadshold(tmp, parent->threadsholdValue->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::gradientFilter()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::gradientEdge(tmp);

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::mopOpen()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::mopOpen(tmp, parent->dilateMask->value(), parent->eroMask->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::mopClose()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::mopClose(tmp, parent->dilateMask->value(), parent->eroMask->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::copyAB()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage original(tmp);

	emit getOriginalImage(original);


	//Filters::subOriginal(tmp, original);

	imageLabel->setPixmap(QPixmap::fromImage(original)); //FIXME

	sendImage(tmp);
}

void ImageViewer::copyBA()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	emit updateImageA(tmp);


	//Filters::subOriginal(tmp, original);

	//	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::saveB()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save File"), QDir::currentPath());
	tmp.save(fileName);
}

void ImageViewer::countObjects()
{
	QImage tmp = imageLabel->pixmap()->toImage();

	int maskValue = parent->objMask->value();
	int cnt = Filters::countObjects(tmp, maskValue);

	QMessageBox::information(this, "Objects Count", "Current objects of mask: " + 
		QString::number(maskValue) + " is: " + QString::number(cnt) );

}

void ImageViewer::gaussianBlur()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	qt_blurImage( tmp, parent->gaussianRadius->value(), true);//, int transposed = 0);
	//Filters::gaussianBlur(tmp, parent->gaussianRadius->value());

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}


void ImageViewer::visWatershed()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage xxx("markers.bmp");
	Filters::visWatershed(tmp, xxx, QString(""));
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}


void ImageViewer::sobel()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::sobelThresh(tmp);
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::canny()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::canny(tmp, parent->cannyRadius->value(), parent->cannyLow->value(), parent->cannyHigh->value());
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::invert()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	tmp.invertPixels();

	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::LoG()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::LoG(tmp, parent->logBlurRadius->value(), parent->logKernelRadius->value());
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::zeroCrossing()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::zeroCrossing(tmp, parent->logBlurRadius->value(), parent->logKernelRadius->value());
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::multiply()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage original(tmp);

	emit getOriginalImage(original);

	Filters::multiply(original, tmp);
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::nullify()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	Filters::nullify(tmp);
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}

void ImageViewer::combine()
{
	QImage tmp = imageLabel->pixmap()->toImage();
	QImage original(tmp);

	emit getOriginalImage(original);

	Filters::combine(original, tmp);
	imageLabel->setPixmap(QPixmap::fromImage(tmp)); //FIXME

	sendImage(tmp);
}


Images::Images(QWidget * parent)
{
	view1 = new ImageViewer(this);
	view2 = new ImageViewer(this);

	hbox1 = new QHBoxLayout(this);

	buttonSNPNoise = new QPushButton("SaltnPeper");
	buttonGaussNoise = new QPushButton("Gauss Noise");
	buttonLPF = new QPushButton("LPF");
	buttonHPF = new QPushButton("HPF");
	buttonBPF = new QPushButton("BPF");
	buttonBSF = new QPushButton("BSF");
	equlise = new QPushButton("Equilise");
	dilate = new QPushButton("Dilate");
	dilateMask = new QSpinBox;
	dilateMask->setValue(3);
	dilateMask->setPrefix("mask: ");
	eroMask = new QSpinBox;
	eroMask->setValue(3);
	eroMask->setPrefix("mask: ");
	erosion = new QPushButton("Erosion");
	dilDifEro = new QPushButton("Dilate-Erosion");
	laplace = new QPushButton("Laplace");
	subBfromA = new QPushButton("Left - Right"); //B-A
	subAfromB = new QPushButton("Right - Left"); //A-B
	gradientFilter = new QPushButton("Gradient");
	morphOpen = new QPushButton("MOP OPEN");
	threadshold = new QPushButton("Threadshold");
	threadsholdValue = new QSpinBox;
	threadsholdValue->setMinimum(0);
	threadsholdValue->setMaximum(255);
	threadsholdValue->setValue(127);
	threadsholdValue->setPrefix("th: ");

	morphClose = new QPushButton("MOP CLOSE");
	spinboxFreq = new QDoubleSpinBox;
	spinboxFreq->setDecimals(4);
	spinboxFreq->setValue(0.1);
	spinboxFreq->setSingleStep(0.01);
	spinboxFreq->setPrefix("fc1: ");
	spinboxFreq->setMaximum(0.5);

	spinboxFreq2 = new QDoubleSpinBox;
	spinboxFreq2->setDecimals(4);
	spinboxFreq2->setValue(0.3);
	spinboxFreq2->setSingleStep(0.01);
	spinboxFreq2->setPrefix("fc2: ");
	spinboxFreq2->setMaximum(0.5);

	snpPerc = new QDoubleSpinBox;
	snpPerc->setDecimals(4);
	snpPerc->setValue(0.03);
	snpPerc->setSingleStep(0.01);
	snpPerc->setPrefix("snp: ");
	snpPerc->setMaximum(1.0);
	snpPerc->setMinimum(0);

	gaussRazb = new QDoubleSpinBox;
	gaussRazb->setValue(10);
	gaussRazb->setPrefix("razb: ");
	gaussRazb->setMinimum(0);

	gaussianRadius = new QDoubleSpinBox;
	gaussianRadius->setValue(50);
	gaussianRadius->setSingleStep(0.01);
	gaussianRadius->setMaximum(9000.0);
	gaussianBlur = new QPushButton("Blur");

	visWater = new QPushButton("Vis Watershed");
	sobelBtn = new QPushButton("Sobel");
	
	logBlurRadius = new QDoubleSpinBox;
	logBlurRadius->setValue(0.5);
	logBlurRadius->setPrefix("sigma: ");
	logKernelRadius = new QSpinBox;
	logKernelRadius->setValue(7);
	logKernelRadius->setPrefix("kernel: ");
	logBtn = new QPushButton("LoG");

	zeroCrossBtn = new QPushButton("Zero crossing");
	multBtn = new QPushButton("Multiply");
	nullifyBtn = new QPushButton("Nullify");
	combineBtn = new QPushButton("Combine");
	equliseBox = new QCheckBox("bg");
	equliseBox->setChecked(false);

	cannyBtn = new QPushButton("Canny");
	cannyRadius = new QDoubleSpinBox;
	cannyRadius->setPrefix("radius: ");
	cannyRadius->setValue(1);
	cannyRadius->setSingleStep(0.1);
	cannyRadius->setMaximum(100);
	cannyLow = new QDoubleSpinBox;
	cannyLow->setPrefix("low: ");
	cannyLow->setValue(1.5);
	cannyLow->setSingleStep(0.1);
	cannyLow->setMaximum(100);
	cannyHigh = new QDoubleSpinBox;
	cannyHigh->setPrefix("high: ");
	cannyHigh->setValue(2);
	cannyHigh->setSingleStep(0.1);
	cannyHigh->setMaximum(100);


	invertBtn = new QPushButton("Invert");

	copyAB = new QPushButton("Copy Left->Right"); //A->B
	copyBA = new QPushButton("Copy Right->Left"); //B->A
	saveB = new QPushButton("Save Right"); //B->A

	objMask = new QSpinBox;
	objMask->setValue(9);
	objMask->setPrefix("mask: ");
	count = new QPushButton("Count"); //B->A


	QVBoxLayout *vbox = new QVBoxLayout;
	QVBoxLayout *vboxNoise = new QVBoxLayout;
	QVBoxLayout *vboxFilters1 = new QVBoxLayout;
	QVBoxLayout *vboxFilters2 = new QVBoxLayout;
	QHBoxLayout *hboxFilters = new QHBoxLayout;

	hbox1->addWidget(view1);
	hbox1->addWidget(view2);
	hbox1->addLayout(vbox);


	QGroupBox *grpNoise = new QGroupBox("Noise");
	QGroupBox *grpFilters = new QGroupBox("Filters");

	vbox->addWidget(grpNoise);
	vbox->addWidget(grpFilters);
	vbox->addStretch(1);

	hboxFilters->addLayout(vboxFilters1);
	hboxFilters->addLayout(vboxFilters2);
	grpFilters->setLayout(hboxFilters);

	vboxFilters1->addWidget(buttonLPF);
	vboxFilters1->addWidget(buttonHPF);
	vboxFilters1->addWidget(buttonBPF);
	vboxFilters1->addWidget(buttonBSF);
	vboxFilters1->addWidget(spinboxFreq);
	vboxFilters1->addWidget(spinboxFreq2);
	vboxFilters1->addWidget(equlise);
	vboxFilters1->addWidget(dilate);
	vboxFilters1->addWidget(dilateMask);
	vboxFilters1->addWidget(erosion);
	vboxFilters1->addWidget(eroMask);
	vboxFilters1->addWidget(dilDifEro);
	vboxFilters1->addWidget(laplace);
	vboxFilters1->addWidget(subBfromA);
	vboxFilters1->addWidget(subAfromB);
	vboxFilters1->addWidget(invertBtn);
	vboxFilters1->addWidget(equliseBox);

	vboxFilters2->addWidget(gradientFilter);
	vboxFilters2->addWidget(threadshold);
	vboxFilters2->addWidget(threadsholdValue);
	vboxFilters2->addWidget(morphOpen);
	vboxFilters2->addWidget(morphClose);
	vboxFilters2->addWidget(copyAB);
	vboxFilters2->addWidget(copyBA);
	vboxFilters2->addWidget(saveB);
	vboxFilters2->addWidget(count);
	vboxFilters2->addWidget(objMask);

	vboxFilters2->addWidget(gaussianBlur);
	vboxFilters2->addWidget(gaussianRadius);
	vboxFilters2->addWidget(visWater);
	vboxFilters2->addWidget(sobelBtn);

	vboxFilters2->addWidget(cannyBtn);
	vboxFilters2->addWidget(cannyRadius);
	vboxFilters2->addWidget(cannyLow);
	vboxFilters2->addWidget(cannyHigh);

	vboxFilters2->addWidget(logBtn);
	vboxFilters2->addWidget(logBlurRadius);
	vboxFilters2->addWidget(logKernelRadius);
	vboxFilters2->addWidget(zeroCrossBtn);
	vboxFilters2->addWidget(multBtn);
	vboxFilters2->addWidget(nullifyBtn);
	vboxFilters2->addWidget(combineBtn);

	grpNoise->setLayout(vboxNoise);
	vboxNoise->addWidget(buttonSNPNoise);
	vboxNoise->addWidget(snpPerc);
	vboxNoise->addWidget(buttonGaussNoise);
	vboxNoise->addWidget(gaussRazb);

	setLayout(hbox1);

}
