#include <QtGui>

#include "noise.h"

//! [0]
Noise::Noise(QWidget *parent)
{

}


void Noise::saltPepper(QImage &img, double percent)
{
	int width = img.width();
	int height = img.height();

	int n = (int)(percent * width * height);
	int xmin = 0;
	int xmax = width;
	int ymin = 0;
	int ymax = height;
	int rx, ry;

	srand(time(NULL));

	for(int i = 0; i < n/2; i++)
        {
			rx = myrandom(xmin, xmax);
			ry = myrandom(ymin, ymax);
			img.setPixel(rx, ry, QColor(255, 255, 255).rgba());

			rx = myrandom(xmin, xmax);
			ry = myrandom(ymin, ymax);
			img.setPixel(rx, ry, QColor(0, 0, 0).rgba());
        }

}

float Noise::gauss (float mu, float sigma2)
{
	double u, v;
	while ((u = rand () / (double) RAND_MAX) == 0.0);
	v = rand () / (double) RAND_MAX;
	/* FIXME: rounding bias here! */
	return mu + sigma2 * sqrt (-2 * log (u)) * cos (2 * M_PI * v);
}

void Noise::gaussNoise(QImage &img, double gaussrazbros)
{
	int width = img.width();
	int height = img.height();

	srand(time(NULL));

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			int r = QColor(img.pixel(i, j)).red();
			r = gauss(r, gaussrazbros);
			if (r > 255)
				r = 255;
			else if (r < 0)
				r = 0;

			img.setPixel(i, j, QColor(r, r, r).rgba());
		}

}
