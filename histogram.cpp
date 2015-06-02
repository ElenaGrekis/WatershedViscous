#include <QtGui>

#include "histogram.h"

//! [0]
Histogram::Histogram(QWidget *parent, QImage *img)
{
	setHistogram = false;
	histogrammArray = new int[256];
}

void Histogram::buildReadHistogramm()
{
	int currsumm = 0;
	for (int i = 0; i < 256; i++)
	{
		currsumm += histogrammArray[i];
		PhistogrammArray[i] = currsumm;
	}

	double divisor = currsumm / 255.0;
	for (int i = 0; i < 256; i++)
	{
		PhistogrammArray[i] = (int)(PhistogrammArray[i] / divisor);
	}
}


void Histogram::readHistogramm(QImage &img)
{
	for(int i = 0; i < 256; i++)
		histogrammArray[i] = 0;

	width = img.width();
	height = img.height();

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			int tmp = QColor(img.pixel(i, j)).red();
			histogrammArray[QColor(img.pixel(i, j)).red()]++;
		}

		setHistogram = true;

		adjustSize();
		update();
}

void Histogram::paintEvent(QPaintEvent *event)
{
	if(setHistogram)
	{

		QPixmap px(width, height);
		px.fill();


		QPainter p2(&px);
		p2.setPen(Qt::black);

		int max = histogrammArray[0];
		for(int i = 0; i < 256; i++)
			if(histogrammArray[i] > max)
				max = histogrammArray[i];

		float DeltaHeight = (float)height / (float)max;
		float DeltaWidth = (float)width / 256;

		for (int i = 0; i < 256; i++)
		{
			p2.drawRect((float)i * DeltaWidth, height - (float)histogrammArray[i] * DeltaHeight, DeltaWidth, height);
		}
		p2.end();

		QPainter p(this);
		p.drawPixmap(0, 0, px);
	}

}

Histograms::Histograms(QWidget *parent)
{
	hbox = new QHBoxLayout(this);

	hbox->addWidget(&hist1);
	hbox->addWidget(&hist2);

	setLayout(hbox);
}