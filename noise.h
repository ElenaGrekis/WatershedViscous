#ifndef NOISE_H
#define NOISE_H

#include <QWidget>
#include <ctime>

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
class Noise : public QWidget
{
    Q_OBJECT

	static inline int myrandom(int min, int max)
	{
		return min + rand() % (max-min);
	}
	
	static inline int compare (const void * a, const void * b)
	{
		return ( *(int*)a - *(int*)b );
	}

	static float gauss (float mu, float sigma2);
public:
    Noise(QWidget *parent = 0);
	static void saltPepper(QImage &img, double percent);
	static void gaussNoise(QImage &img, double gaussrazbros);


};
//! [0]
#endif
