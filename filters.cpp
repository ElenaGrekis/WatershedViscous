#include <QtGui>

#include "filters.h"
#include "laplacianZC.h"
#include "qtconvert.h"
//! [0]
Filters::Filters(QWidget *parent)
{

}

#define MAX_BRIGHTNESS  255

void Filters::main_filter(QImage &img, Filters2 fil, double fc, double fc2)
{
	int width = img.width();
	int height = img.height();


	double filter_dt = 1.0; //?
	double *arr3 = nullptr;
	if(fil == LPF)
		arr3 = lp(fc, m, filter_dt);
	else if(fil == HPF)
		arr3 = hp(fc, m, filter_dt);
	else if(fil == BPF)
		arr3 = bp(fc, fc2, m, filter_dt);
	else if(fil == BSF)
		arr3 = bs(fc, fc2, m, filter_dt);

	for(int i = 0; i < width; i++)
	{
		double *arr2 = new double[height];
		for(int j = 0; j < height; j++)
		{
			arr2[j] = QColor(img.pixel(i, j)).red();
		}

		double *arr4 = svertka(arr2, height, arr3, 2*m + 1);

		for(int j = 0; j < height; j++)
		{
			/*if(arr4[j] > 255)
			arr4[j] = 255;
			if(arr4[j] < 0)
			arr4[j] = 255;*/
			img.setPixel(i, j, QColor((int)arr4[j], (int)arr4[j], (int)arr4[j]).rgba());
		}

		delete [] arr2;
		//delete [] arr4;
	}


	for(int j = 0; j < height; j++)
	{
		double *arr2 = new double[width];
		for(int i = 0; i < width; i++)
		{
			arr2[i] = QColor(img.pixel(i, j)).red();
		}

		double *arr4 = svertka(arr2, width, arr3, 2*m + 1);

		for(int i = 0; i < width; i++)
		{
			img.setPixel(i, j, QColor((int)arr4[i], (int)arr4[i], (int)arr4[i]).rgba());
		}

		delete [] arr2;
		//delete [] arr4;
	}



	delete [] arr3;

}

double *Filters::svertka(double *val, int val_size, double *filt, int filt_size)
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
		return &arr_tmp[L/2];
}


double * Filters::lp(double f0, const int M, double dt) 
{
	double d[] = {0.35577019, 0.2436983, 0.07211497, 0.00630165 };
	double *wlp = new double [2 * m + 1];
	double tmp[m + 1];
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
	for (int i = 0; i <= m; i++)
		tmp[i] /= sumg;


	for (int i = 0; i < m; i++)
		wlp[m - 1 - i] = wlp[m + 1 + i] = tmp[i + 1];
	wlp[m] = tmp[0];
	return wlp;


}

double *Filters::hp(double f0, const int M, double dt)
{
	double *whp = new double [2 * m + 1];
	double *wlp = lp(f0, m, dt);
	for (int i = 1; i < m + 1; i++)
		whp[m + i] = whp[m - i] = - wlp[m + i];
	whp[m] = 1 - wlp[m];

	delete [] wlp;
	return whp;

}


double* Filters::bp(double f1, double f2, const int M, double dt)
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

double* Filters::bs(double f1, double f2, const int M, double dt)
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



void Filters::equilisation(QImage &img, bool ignoreBackground)
{
	
	int histogrammArray[256] = {0};
	int PhistogrammArray[256] = {0};
	int width = img.width();
	int height = img.height();

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			histogrammArray[QColor(img.pixel(i, j)).red()]++;
		}

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

		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				int r = QColor(img.pixel(i, j)).red();
				int nr = PhistogrammArray[r];
				if(ignoreBackground)
				{
					if(r == 255)
					{}
					else
						img.setPixel(i, j, QColor(nr, nr, nr).rgba());
				}
				else
					img.setPixel(i, j, QColor(nr, nr, nr).rgba());
			}
}


void Filters::fillOne(int **matr, int w, int h)
{
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			matr[i][j] = 1;

}

void Filters::dilate(QImage &img, int mask_size)
{
	int width = img.width();
	int height = img.height();
	int max = 0;

	int **mask = new int*[mask_size];
	for(int i = 0; i < mask_size; i++)
		mask[i] = new int[mask_size];

	fillOne(mask, mask_size, mask_size);


	QImage tmp(img);

	int mask_c = mask_size / 2;

	for(int i = mask_c; i < width - mask_c; i++)
	{
		for(int j = mask_c; j < height - mask_c; j++)
		{
			max = 0;
			for(int a = -mask_c; a <= mask_c; a++)
			{
				for(int b = -mask_c; b <= mask_c; b++)
				{
					if((a+mask_c < mask_size) && (b+mask_c < mask_size) && mask[a+mask_c][b+mask_c] == 1)
					{
						if(QColor(img.pixel(i+a, j+b)).red() > max)
							max = QColor(img.pixel(i+a, j+b)).red();
					}
				}
			}

			tmp.setPixel(i, j, QColor(max, max, max).rgba());
		}
	}

	delete [] mask;

	img = tmp;
}

void Filters::erosion(QImage &img, int mask_size)
{

	int width = img.width();
	int height = img.height();
	int min = 0;


	int **mask = new int*[mask_size];
	for(int i = 0; i < mask_size; i++)
		mask[i] = new int[mask_size];

	fillOne(mask, mask_size, mask_size);


	QImage tmp(img);

	int mask_c = mask_size / 2;

	for(int i = mask_c; i < width - mask_c; i++)
	{
		for(int j = mask_c; j < height - mask_c; j++)
		{
			min = 255;
			for(int a = -mask_c; a <= mask_c; a++)
			{
				for(int b = -mask_c; b <= mask_c; b++)
				{
					if((a+mask_c < mask_size) && (b+mask_c < mask_size) && mask[a+mask_c][b+mask_c] == 1)
					{
						if(QColor(img.pixel(i+a, j+b)).red() < min)
							min = QColor(img.pixel(i+a, j+b)).red();
					}
				}
			}

			tmp.setPixel(i, j, QColor(min, min, min).rgba());
		}
	}

	img = tmp;
}


void Filters::DilDifEro(QImage &img, int dilate_mask, int ero_mask)
{

	int width = img.width();
	int height = img.height();


	QImage dil(img);
	QImage ero(img);

	dilate(dil, dilate_mask);
	erosion(ero, ero_mask);

	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			int res = QColor(dil.pixel(i, j)).red() - QColor(ero.pixel(i, j)).red();
			img.setPixel(i, j, QColor(res, res, res).rgba());
		}

}

void Filters::laplace(QImage &img)
{

	int width = img.width();
	int height = img.height();

	QImage tmp(img);
	tmp.fill(0);

	int weight[3][3] = {{ -1,  -1,  -1 },
	{ -1,  8,  -1 },
	{ -1,  -1,  -1 }};

	double pixel_value;
	double min, max;
	int x, y, i, j; 


	min = DBL_MAX;
	max = -DBL_MAX;
	for (y = 1; y < width - 1; y++) {
		for (x = 1; x < height - 1; x++) {
			pixel_value = 0.0;
			for (j = - 1; j < 2; j++) {
				for (i = -1; i < 2; i++) {
					pixel_value += weight[j + 1][i + 1] * QColor(img.pixel(y + j, x + i)).red();
				}
			}
			if (pixel_value < min) min = pixel_value;
			if (pixel_value > max) max = pixel_value;
		}
	}

	for (y = 1; y < width - 1; y++) {
		for (x = 1; x < height - 1; x++) {
			pixel_value = 0.0;
			for (j = - 1; j < 2; j++) {
				for (i = -1; i < 2; i++) {
					pixel_value += weight[j + 1][i + 1] * QColor(img.pixel(y + j, x + i)).red();
				}
			}

			//pixel_value = MAX_BRIGHTNESS * (pixel_value - min) / (max - min); //!!!!For cool gray
			if (pixel_value > 255)
			{
				pixel_value = 255; 
			}
			else if (pixel_value < 0)
			{ 
				pixel_value = 0; 
			}
			unsigned char px = (unsigned char)pixel_value;
			tmp.setPixel(y, x, QColor(px, px, px).rgba()); 
		}
	}

	img = tmp;
}		


void Filters::subOriginal(QImage &img1, QImage &orig)
{
	int width = img1.width();
	int height = img1.height();


	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			int res = QColor(orig.pixel(i, j)).red() - QColor(img1.pixel(i, j)).red();
			img1.setPixel(i, j, QColor(res, res, res).rgba());
		}

}

int Filters::countObjects(QImage &img, int mask_size)
{

	int width = img.width();
	int height = img.height();
	int min = 0;


	int **mask = new int*[mask_size];
	for(int i = 0; i < mask_size; i++)
		mask[i] = new int[mask_size];

	fillOne(mask, mask_size, mask_size);


	QImage tmp(img);

	int mask_c = mask_size / 2;

	int countObjects = 0;

	for(int i = mask_c; i < width - mask_c; i++)
	{
		for(int j = mask_c; j < height - mask_c; j++)
		{
			min = 255;
			bool same = false;
			for(int a = -mask_c; a <= mask_c; a++)
			{
				for(int b = -mask_c; b <= mask_c; b++)
				{
					if((a+mask_c < mask_size) && (b+mask_c < mask_size) && mask[a+mask_c][b+mask_c] == 1)
					{
						if(QColor(img.pixel(i+a, j+b)).red() == min)
						{
							same = true;
						}
						else
						{
							same = false;
							break;
						}
						
					}
				}
				if(!same)
					break;
			}
			if(same)
				countObjects++;
		}
	}

	return countObjects;
}



int Filters::gradient (const int *data)
{
	const int v_kernel[9] = { 0,  0,  0,
		0,  4, -4,
		0,  0,  0 };
	const int h_kernel[9] = { 0,  0,  0,
		0, -4,  0,
		0,  4,  0 };

	int i;
	int v_grad, h_grad;

	for (i = 0, v_grad = 0, h_grad = 0; i < 9; i++)
	{
		v_grad += v_kernel[i] * data[i];
		h_grad += h_kernel[i] * data[i];
	}

	return  sqrt ((double)v_grad * v_grad  +
		(double)h_grad * h_grad );
}


void Filters::gradientEdge(QImage &img)
{

	int width = img.width();
	int height = img.height();

	QImage tmp(img);
	tmp.fill(0);

	int weight[3][3] = { 0 };

	double pixel_value;
	double min, max;
	int x, y, i, j; 


	min = DBL_MAX;
	max = -DBL_MAX;
	for (y = 1; y < width - 1; y++)
	{
		for (x = 1; x < height - 1; x++) 
		{
			pixel_value = 0.0;
			for (j = - 1; j < 2; j++) 
			{
				for (i = -1; i < 2; i++) 
				{
					weight[j + 1][i + 1] = QColor(img.pixel(y + j, x + i)).red();
				}
			}
			pixel_value = gradient(weight[0]);

			if (pixel_value < min) min = pixel_value;
			if (pixel_value > max) max = pixel_value;
		}
	}

	for (y = 1; y < width - 1; y++) {
		for (x = 1; x < height - 1; x++) {
			pixel_value = 0.0;
			for (j = - 1; j < 2; j++) {
				for (i = -1; i < 2; i++) {
					weight[j + 1][i + 1] = QColor(img.pixel(y + j, x + i)).red();
				}
			}
			pixel_value = gradient(weight[0]);

			if (pixel_value > 255)
			{
				pixel_value = 255; 
			}
			else if (pixel_value < 0)
			{ 
				pixel_value = 0; 
			}
			unsigned char px = (unsigned char)pixel_value;
			tmp.setPixel(y, x, QColor(px, px, px).rgba()); 
		}
	}

	img = tmp;
}		

void Filters::threadshold(QImage &img, int th)
{

	int width = img.width();
	int height = img.height();

	QImage tmp(img);

	int pixel_value = 0;
	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			pixel_value = QColor(tmp.pixel(i, j)).red();

			if(pixel_value < th)
			{
				tmp.setPixel(i, j, QColor(0, 0, 0).rgba());
			}
			else if(pixel_value >= th)
			{
				tmp.setPixel(i, j, QColor(255, 255, 255).rgba());
			}
		}

		img = tmp;
}		


void Filters::mopOpen(QImage &img, int dilate_mask, int ero_mask)
{

	int width = img.width();
	int height = img.height();

	QImage tmp(img);

	erosion(tmp, ero_mask);
	dilate(tmp, dilate_mask);

	img = tmp;
}		


void Filters::mopClose(QImage &img, int dilate_mask, int ero_mask)
{

	int width = img.width();
	int height = img.height();

	QImage tmp(img);

	dilate(tmp, dilate_mask);
	erosion(tmp, ero_mask);

	img = tmp;
}

static void minMax(const QImage &img, int &min, int &max)
{

	int width = img.width();
	int height = img.height();


	min = UCHAR_MAX;
	max = -UCHAR_MAX;


	int pixel_value = 0;
	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			pixel_value = QColor(img.pixel(i, j)).red();

			if (pixel_value < min) 
				min = pixel_value;
			if (pixel_value > max) 
				max = pixel_value;
		}


}

void Filters::gaussianBlur(QImage &img, float radius)
{
//	qt_blurImage( img, radius, true);//, int transposed = 0);
	
	int mWidth = img.width();
	int mHeight = img.height();

	QImage tmp(img);

	float sigma2=radius*radius;
	int size=5; //размер фильтра


	float pixel = 0;
	float sum = 0;

	//blur изображения по x 
	for(int y=0;y<mHeight;y++)
	{
		for(int x=0;x<mWidth;x++)
		{
			sum=0;
			pixel = 0;

			//сложение пикселей
			for(int i=qMax(0,x-size);i<=qMin(mWidth-1,x+size);i++)
			{
				float factor=exp(-(i-x)*(i-x)/(2.0*sigma2));
				sum+=factor;
				pixel += factor * QColor(tmp.pixel(i, y)).red();
			}

			//запись значения пикселя в картинку
				float temp = pixel/sum;
				tmp.setPixel(x, y, QColor(temp, temp, temp).rgba());
		}
	}

	//blur изображения по y
	for(int y=0;y<mHeight;y++)
	{
		for(int x=0;x<mWidth;x++)
		{
			sum = 0;
			pixel = 0;

			//сложение пикселей
			for(int i=qMax(0,y-size);i<=qMin(mHeight-1,y+size);i++)
			{
				float factor=exp(-(i-y)*(i-y)/(2.0*sigma2));
				sum+=factor;

					pixel += factor * QColor(tmp.pixel(x, i)).red();

			}

			//запись значения пикселя в картинку
				float temp = pixel / sum;
				img.setPixel(x, y, QColor(temp, temp, temp).rgba());

		}
	}
}

CvWSNode* icvAllocWSNodes( CvMemStorage* storage )
{
	CvWSNode* n = 0;	
	int i, count = (storage->block_size - sizeof(CvMemBlock))/sizeof(*n) - 1;
	n = (CvWSNode*)cvMemStorageAlloc( storage, count*sizeof(*n) );
	for( i = 0; i < count-1; i++ )
		n[i].next = n + i + 1;
	n[count-1].next = 0;
	return n;
}

//маркерный водораздел с вязкостью
void Filters::visWatershed(QImage &imgTemp, QImage &markersTemp, QString &path)
{
	int sz = 2; 
	double C = 0.02;
	//путь к картинке
	QString picturePath = QFileDialog::getOpenFileName(0,
		tr("Open Picture"), QDir::currentPath());

	//путь к маркеру маркера
	QString maskPath = QFileDialog::getOpenFileName(0,
		tr("Open Mask"), QDir::currentPath());
	if (maskPath.isEmpty()) 
		return;

	//загрузка маркера
	IplImage *i1 = cvCreateImage(cvSize(imgTemp.width(), imgTemp.height()), IPL_DEPTH_8U, 1);
	i1 = cvLoadImage(maskPath.toLocal8Bit().constData(), 0);
	IplImage *i2 = cvCreateImage(cvGetSize(i1), IPL_DEPTH_32S, 1);
	cvConvertScale(i1, i2);
	cv::Mat mat2(i2);
	CvMat mmm2 = mat2;

	//повторная загрузка картинки
	IplImage *i0 = cvLoadImage(picturePath.toLocal8Bit().constData(), 1);
	cv::Mat mat(i0);

	CvMat mmm = mat;


	IplImage *clone = cvLoadImage(picturePath.toLocal8Bit().constData(), 1);
	cvZero(clone);

	const int IN_QUEUE = -2;
	const int WSHED = -1;
	const int NQ = 256;

	const int thres = 256;
	const int block = (2*sz+1)*(2*sz+1);
	CvMemStorage* storage = 0;

	CvMat sstub, *src;
	CvMat dstub, *dst;
	CvSize size;
	CvWSNode* free_node = 0, *node;
	CvWSQueue q[NQ];
	int active_queue;
	int i, j;
	int db, dg, dr;
	int* mask;
	uchar* img;
	int mstep, istep;
	int subs_tab[513];


#define ws_max(a,b) ((b) + subs_tab[(a)-(b)+NQ])
#define ws_min(a,b) ((a) - subs_tab[(a)-(b)+NQ])
#define ws_push(idx,mofs,iofs)  \
	{                               \
	if( !free_node )            \
	free_node = icvAllocWSNodes( storage );\
	node = free_node;           \
	free_node = free_node->next;\
	node->next = 0;             \
	node->mask_ofs = mofs;      \
	node->img_ofs = iofs;       \
	if( q[idx].last )           \
	q[idx].last->next=node; \
	else                        \
	q[idx].first = node;    \
	q[idx].last = node;         \
	}

#define ws_pop(idx,mofs,iofs)   \
	{                               \
	node = q[idx].first;        \
	q[idx].first = node->next;  \
	if( !node->next )           \
	q[idx].last = 0;        \
	node->next = free_node;     \
	free_node = node;           \
	mofs = node->mask_ofs;      \
	iofs = node->img_ofs;       \
	}

#define c_diff(ptr1,ptr2,diff)      \
	{                                   \
	db = abs((ptr1)[0] - (ptr2)[0]);\
	dg = abs((ptr1)[1] - (ptr2)[1]);\
	dr = abs((ptr1)[2] - (ptr2)[2]);\
	diff = ws_max(db,dg);           \
	diff = ws_max(diff,dr);         \
	assert( 0 <= diff && diff <= 255 ); \
	}


	src = &mmm;
	dst = &mmm2;

	size = cvGetMatSize(src);

	storage = cvCreateMemStorage();

	istep = src->step;
	img = src->data.ptr;
	mstep = dst->step / sizeof(mask[0]);
	mask = dst->data.i;

	memset( q, 0, NQ*sizeof(q[0]) );

	for( i = 0; i < 256; i++ )
		subs_tab[i] = 0;
	for( i = 256; i <= 512; i++ )
		subs_tab[i] = i - 256;

	//рисование рамки из пикселей вокруг изображения (border of dummy "watershed")
	//записывается число -1 (255). WSHED - константа.
	for( i = 0; i < sz; i++ ) 
		for( j = 0; j < size.width; j++ ){
			mask[j+i*mstep] = mask[j + mstep*(size.height-i-1)] = WSHED;
		}

		mask += (sz-1)*mstep;
		img += (sz-1)*istep;

		//начальная фаза. все соседние пиксели каждого маркера помещаются в упорядоченную очередь
		//определяются начальные границы бассейнов
		for( i = sz; i < size.height-sz; i++ )
		{
			img += istep; mask += mstep;
			for( int n=0; n < sz; n++ )
				mask[0+n] = mask[size.width-n-1] = WSHED;

			for( j = sz; j < size.width-sz; j++ )
			{
				int* m = mask + j;
				if( m[0] < 0 ) m[0] = 0;
				if( m[0] == 0 && (m[-1] > 0 || m[1] > 0 || m[-mstep] > 0 || m[mstep] > 0) )
				{
					uchar* ptr = img + j*3;
					int idx = 256, t;
					if( m[-1] > 0 )
						c_diff( ptr, ptr - 3, idx );
					if( m[1] > 0 )
					{
						c_diff( ptr, ptr + 3, t );
						idx = ws_min( idx, t );
					}
					if( m[-mstep] > 0 )
					{
						c_diff( ptr, ptr - istep, t );
						idx = ws_min( idx, t );
					}
					if( m[mstep] > 0 )
					{
						c_diff( ptr, ptr + istep, t );
						idx = ws_min( idx, t );
					}
					assert( 0 <= idx && idx <= 255 );
					ws_push( idx, i*mstep + j, i*istep + j*3 );
					m[0] = IN_QUEUE;
				}
			}
		}

		//поиск первой непустой очереди
		for( i = 0; i < NQ; i++ )
			if( q[i].first )
				break;

		//если нет маркеров, то выход, т.е. :
		//число очередей - 255, если NQ будет 256, то выход
		if( i == NQ )
			exit;

		active_queue = i;
		img = src->data.ptr;
		mask = dst->data.i;

		//рекурсивное заполнение бассейнов
		for(;;)
		{
			int mofs, iofs;
			int lab = 0, t, t1, b;
			int* m;
			uchar* ptr;

			if( q[active_queue].first == 0 )
			{
				for( i = active_queue+1; i < NQ; i++ )
					if( q[i].first )
						break;

				if( i == NQ )
					break;
				active_queue = i;
			}

			ws_pop( active_queue, mofs, iofs );

			m = mask + mofs;
			ptr = img + iofs;
			t = m[-1];
			if( t > 0 ) lab = t;
			t = m[1];
			if( t > 0 )
				if( lab == 0 ) lab = t;
				else if( t != lab ) lab = WSHED;
				t = m[-mstep];
				if( t > 0 )
					if( lab == 0 ) lab = t;
					else if( t != lab ) lab = WSHED;
					t = m[mstep];
					if( t > 0 )
						if( lab == 0 ) lab = t;
						else if( t != lab ) lab = WSHED;

						assert( lab != 0 );
						m[0] = lab;
						if( lab == WSHED )
							continue;

						if( m[-1] == 0 )
						{
							c_diff( ptr, ptr - 3, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr-3, ptr-i*3+j*istep, t1);
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);   
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs - 1, iofs - 3 );
							active_queue = ws_min( active_queue, t );
							m[-1] = IN_QUEUE;
						}
						if( m[1] == 0 )
						{
							c_diff( ptr, ptr + 3, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr+3, ptr+i*3+j*istep, t1);
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs + 1, iofs + 3 );
							active_queue = ws_min( active_queue, t );
							m[1] = IN_QUEUE;
						}
						if( m[-mstep] == 0 )
						{
							c_diff( ptr, ptr-istep, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr-istep, ptr+j*3-i*istep, t1 );
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs - mstep, iofs - istep );
							active_queue = ws_min( active_queue, t );
							m[-mstep] = IN_QUEUE;
						}
						if( m[mstep] == 0 )
						{
							c_diff( ptr, ptr+istep, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr+istep, ptr+j*3+i*istep, t1 );
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs + mstep, iofs + istep );
							active_queue = ws_min( active_queue, t );
							m[mstep] = IN_QUEUE;
						}
		}
		cvReleaseMemStorage( &storage );

		//получение результирующего изображения 
		//обработанного водоразделом
		for( int i = 0; i < i2->height; i++ )
			for( int j = 0; j < i2->width; j++ )
			{
				int idx = CV_IMAGE_ELEM( i2, int, i, j ); //i2 - финальное изображение с бассейнами водоразделов
				uchar* dst1 = &CV_IMAGE_ELEM( clone, uchar, i, j*3 );

				if( idx == -1 )
				{
					dst1[0] = (uchar)0; dst1[1] = (uchar)255; dst1 [2] = (uchar)0;
					if (i+1< i2->height && i-1>=0 && j+1<i2->width && j-1>=0) 
					{
						for (int n=-1; n<=1; n++)
							for (int m=-1; m<=1; m++)
							{
								uchar* dst2 = &CV_IMAGE_ELEM( clone, uchar, i+n, (j+m)*3 );
								dst2[0] = (uchar)0; dst2[1] = (uchar)0; dst2[2] = (uchar)255;
								//отмечаем красную рамочку на исходном изображении
								imgTemp.setPixel(j+m, i+n, QColor(255, 0, 0).rgba());
							}
					}
				}
			}
			/*
		for( int i = 10; i < i2->height-10; i++ )
			for( int j = 10; j < i2->width-10;  )
			{
				if(QColor(255, 0, 0) == imgTemp.pixel(j, i))
				{
					for( ; j < i2->width-10; j++ )
					{
						imgTemp.setPixel(j, i, QColor(0, 0, 0).rgba());
						if(QColor(255, 0, 0) == imgTemp.pixel(j, i))
							break;
					}
				}
				else
					j++;
			}
		*/
		cvConvertScale(i2, i1, 1, 0);
		cvSaveImage("final_markers.png", i1);

}

int** setImage(int numRows, int numCols)
{
	int** imgTemp = new int*[numRows];

	for(unsigned int i = 0; i < numRows; ++i)
	{
		imgTemp[i] = new int [numCols];
	}

	return imgTemp;
}

void delImage(int **image, int numRows)
{
	for(unsigned int i = 0; i < numRows; ++i)
	{
		delete [] image[i];
	}

	delete [] image;
}

//собель
void Filters::sobelThresh(QImage &pic)
{

	int numRows = pic.height(),numCols = pic.width();
   
	char ci;

	//Маски для свертки 3x3 по направлениям x и y 
	int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
	int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};

	int maxival;
	int i,j,p,q,mr,sum1,sum2;

	int **tempx = setImage(numRows, numCols);
	int **tempy = setImage(numRows, numCols);
	int **grad = setImage(numRows, numCols);

	mr = 1;

	for (i=mr;i<numRows-mr;i++)
	{ 
		for (j=mr;j<numCols-mr;j++)
		{
			sum1 = 0;
			sum2 = 0;
			//Свертка на регионах 3x3
			for (p=-mr;p<=mr;p++)
			{
				for (q=-mr;q<=mr;q++)
				{
					sum1 += QColor(pic.pixel(j+q, i+p)).red() * maskx[p+mr][q+mr];
					sum2 += QColor(pic.pixel(j+q, i+p)).red() * masky[p+mr][q+mr];
				}
			}
			//Сохранение результатов свертки во временных изображениях
			tempx[i][j] = sum1;
			tempy[i][j] = sum2;
		}
	}

	//Вычисление порога на градиентном изображении используя массивы tempx tempy
	maxival = 0;
	for (i=mr;i<numRows-mr;i++)
	{ 
		for (j=mr;j<numCols-mr;j++)
		{
			//Вычисление сглаженного изображения
			grad[i][j] = (int) (sqrt((double)((tempx[i][j]*tempx[i][j]) 
				+ (tempy[i][j]*tempy[i][j]))));

			if (grad[i][j] > maxival)
				maxival = grad[i][j];

		}
	}

	// Сохранение результата работы метода в оригинальном изображении pic
	for (i=0;i<numRows;i++)
	{ 
		for (j=0;j<numCols;j++)
		{
			pic.setPixel(j, i, QColor(grad[i][j], grad[i][j], grad[i][j]).rgba());
		}
	}

	//Удаление динамической памяти
	delImage(tempx, numRows);
	delImage(tempy, numRows);
	delImage(grad, numRows);
}

//Canny создание и инициализация памяти под статические переменные
float Filters::GAUSSIAN_CUT_OFF = 0.005f;
//Множитель для нижнего и верхнего порога
float Filters::MAGNITUDE_SCALE = 100.0f;
//Максимальная величина порога
float Filters::MAGNITUDE_LIMIT = 1000.0f;
int Filters::MAGNITUDE_MAX = (int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT);
int	Filters::gaussianKernelWidth = 16;

int Filters::width;
int Filters::height;
int Filters::picsize;
int *Filters::data;
int *Filters::magnitude;

float *Filters::xConv;
float *Filters::yConv;
float *Filters::xGradient;
float *Filters::yGradient;

//radius - размеря ядра для гауссового фильтра, low и high - верхняя и нижняя границы
void Filters::canny(QImage &img, float radius, float low, float high)
{
	width = img.width();
	height  = img.height();
	picsize = width * height;


	data = new int[picsize];
	magnitude = new int[picsize];

	xConv = new float[picsize];
	yConv = new float[picsize];
	xGradient = new float[picsize];
	yGradient = new float[picsize];

	//Преобразование картинки в массив байт для работы фильтра
	unsigned int byte_count = 0;
	for(int i = 0; i < height; i++)
	{
		for(int j = 0; j < width; j++)
		{
			data[byte_count] = QColor(img.pixel(j, i)).red();
			byte_count++;
		}
	}
	
	//Сглаживание и поиск градиентов
	computeGradients(radius, gaussianKernelWidth);
	int lowP = round(low * MAGNITUDE_SCALE);
	int highP = round( high * MAGNITUDE_SCALE);
	performHysteresis(lowP, highP);
	thresholdEdges();

	delete [] yGradient;
	delete [] xGradient;
	delete [] yConv;
	delete [] xConv;
	delete [] magnitude;

	byte_count = 0;
	for(int i = 0; i < height; i++)
	{
		for(int j = 0; j < width; j++)
		{
			img.setPixel(j, i, QColor(data[byte_count], data[byte_count], data[byte_count]).rgba());
			byte_count++;
		}
	}
	delete [] data;

	img.invertPixels(); //Инвертирование картинки, чтобы были контуры - черные, а фон - белый
}


void Filters::computeGradients(float kernelRadius, int kernelWidth) {

	//маски для свертки гауссом
	float *kernel = new float[kernelWidth];
	float *diffKernel = new float[kernelWidth];
	int kwidth;
	for (kwidth = 0; kwidth < kernelWidth; kwidth++) 
	{
		float g1 = gaussian(kwidth, kernelRadius);
		if (g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2) break;
		float g2 = gaussian(kwidth - 0.5f, kernelRadius);
		float g3 = gaussian(kwidth + 0.5f, kernelRadius);
		kernel[kwidth] = (g1 + g2 + g3) / 3.0f / (2.0f * (float) M_PI * kernelRadius * kernelRadius);
		diffKernel[kwidth] = g3 - g2;
	}

	int initX = kwidth - 1;
	int maxX = width - (kwidth - 1);
	int initY = width * (kwidth - 1);
	int maxY = width * (height - (kwidth - 1));

	//свертка по двум направлениям x и y
	for (int x = initX; x < maxX; x++) 
	{
		for (int y = initY; y < maxY; y += width) 
		{
			int index = x + y;
			float sumX = data[index] * kernel[0];
			float sumY = sumX;
			int xOffset = 1;
			int yOffset = width;
			for(; xOffset < kwidth ;) {
				sumY += kernel[xOffset] * (data[index - yOffset] + data[index + yOffset]);
				sumX += kernel[xOffset] * (data[index - xOffset] + data[index + xOffset]);
				yOffset += width;
				xOffset++;
			}

			yConv[index] = sumY;
			xConv[index] = sumX;
		}

	}

	for (int x = initX; x < maxX; x++)
	{
		for (int y = initY; y < maxY; y += width)
		{
			float sum = 0.0f;
			int index = x + y;
			for (int i = 1; i < kwidth; i++)
				sum += diffKernel[i] * (yConv[index - i] - yConv[index + i]);

			xGradient[index] = sum;
		}

	}

	for (int x = kwidth; x < width - kwidth; x++) {
		for (int y = initY; y < maxY; y += width) {
			float sum = 0.0f;
			int index = x + y;
			int yOffset = width;
			for (int i = 1; i < kwidth; i++) {
				sum += diffKernel[i] * (xConv[index - yOffset] - xConv[index + yOffset]);
				yOffset += width;
			}

			yGradient[index] = sum;
		}

	}


	initX = kwidth;
	maxX = width - kwidth;
	initY = width * kwidth;
	maxY = width * (height - kwidth);
	for (int x = initX; x < maxX; x++) {
		for (int y = initY; y < maxY; y += width) {
			int index = x + y;
			int indexN = index - width;
			int indexS = index + width;
			int indexW = index - 1;
			int indexE = index + 1;
			int indexNW = indexN - 1;
			int indexNE = indexN + 1;
			int indexSW = indexS - 1;
			int indexSE = indexS + 1;

			float xGrad = xGradient[index];
			float yGrad = yGradient[index];
			float gradMag = hypot(xGrad, yGrad);

			//Выполнение подавления немаксимумов.
			//hypot - эвклидово расстояние.
			//sqrt(x*x + y*y)
			float nMag = hypot(xGradient[indexN], yGradient[indexN]);
			float sMag = hypot(xGradient[indexS], yGradient[indexS]);
			float wMag = hypot(xGradient[indexW], yGradient[indexW]);
			float eMag = hypot(xGradient[indexE], yGradient[indexE]);
			float neMag = hypot(xGradient[indexNE], yGradient[indexNE]);
			float seMag = hypot(xGradient[indexSE], yGradient[indexSE]);
			float swMag = hypot(xGradient[indexSW], yGradient[indexSW]);
			float nwMag = hypot(xGradient[indexNW], yGradient[indexNW]);
			float tmp;
			/*
			* An explanation of what's happening here, for those who want
			* to understand the source: This performs the "non-maximal
			* supression" phase of the Canny edge detection in which we
			* need to compare the gradient magnitude to that in the
			* direction of the gradient; only if the value is a local
			* maximum do we consider the point as an edge candidate.
			* 
			* We need to break the comparison into a number of different
			* cases depending on the gradient direction so that the
			* appropriate values can be used. To avoid computing the
			* gradient direction, we use two simple comparisons: first we
			* check that the partial derivatives have the same sign (1)
			* and then we check which is larger (2). As a consequence, we
			* have reduced the problem to one of four identical cases that
			* each test the central gradient magnitude against the values at
			* two points with 'identical support'; what this means is that
			* the geometry required to accurately interpolate the magnitude
			* of gradient function at those points has an identical
			* geometry (upto right-angled-rotation/reflection).
			* 
			* When comparing the central gradient to the two interpolated
			* values, we avoid performing any divisions by multiplying both
			* sides of each inequality by the greater of the two partial
			* derivatives. The common comparand is stored in a temporary
			* variable (3) and reused in the mirror case (4).
			* 
			*/
			if (xGrad * yGrad <= (float) 0 /*(1)*/
				? abs(xGrad) >= abs(yGrad) /*(2)*/
				? (tmp = abs(xGrad * gradMag)) >= abs(yGrad * neMag - (xGrad + yGrad) * eMag) /*(3)*/
				&& tmp > abs(yGrad * swMag - (xGrad + yGrad) * wMag) /*(4)*/
				: (tmp = abs(yGrad * gradMag)) >= abs(xGrad * neMag - (yGrad + xGrad) * nMag) /*(3)*/
				&& tmp > abs(xGrad * swMag - (yGrad + xGrad) * sMag) /*(4)*/
				: abs(xGrad) >= abs(yGrad) /*(2)*/
				? (tmp = abs(xGrad * gradMag)) >= abs(yGrad * seMag + (xGrad - yGrad) * eMag) /*(3)*/
				&& tmp > abs(yGrad * nwMag + (xGrad - yGrad) * wMag) /*(4)*/
				: (tmp = abs(yGrad * gradMag)) >= abs(xGrad * seMag + (yGrad - xGrad) * sMag) /*(3)*/
				&& tmp > abs(xGrad * nwMag + (yGrad - xGrad) * nMag) /*(4)*/
				) 
			{
					magnitude[index] = gradMag >= MAGNITUDE_LIMIT ? MAGNITUDE_MAX : (int) (MAGNITUDE_SCALE * gradMag);

			} 
			else 
			{
				magnitude[index] = 0;
			}
		}
	}

	delete 	[] kernel;
	delete 	[] diffKernel;
}


//Выполнение двойной пороговой фильтрации по границам low и high
void Filters::performHysteresis(int low, int high) 
{
	for(unsigned int i = 0; i < picsize; i++)
		data[i] = 0;

	int offset = 0;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++) 
		{
			if (data[offset] == 0 && magnitude[offset] >= high) 
			{
				//трассировка областей неоднозначности
				follow(x, y, offset, low);
			}
			offset++;
		}
	}
}

//трассировка областей неоднозначности
void Filters::follow(int x1, int y1, int i1, int threshold) 
{
	int x0 = x1 == 0 ? x1 : x1 - 1;
	int x2 = x1 == width - 1 ? x1 : x1 + 1;
	int y0 = y1 == 0 ? y1 : y1 - 1;
	int y2 = y1 == height -1 ? y1 : y1 + 1;

	data[i1] = magnitude[i1];
	for (int x = x0; x <= x2; x++) {
		for (int y = y0; y <= y2; y++) {
			int i2 = x + y * width;
			if ((y != y1 || x != x1)
				&& data[i2] == 0 
				&& magnitude[i2] >= threshold) {
					follow(x, y, i2, threshold);
					return;
			}
		}
	}
}

//Бинаризация изображения
void Filters::thresholdEdges()
{
	for (int i = 0; i < picsize; i++) 
	{
		data[i] = data[i] > 0 ? 255 : 0;
	}
}
void Filters::LoG(QImage &img, float blurRadius, float kernelRadius)
{
	//Mat src = QImageToCvMat(img);;//imread("C:\\qwt-6.1.0\\examples\\imageviewer4\\05-11-2014_13-13-24\\a12-103h_original.jpg", 0);//QImageToCvMat(img).clone();
	Mat src = image2Mat(img, CV_8UC1, getColorOrderOfRGB32Format());
	GaussianBlur( src, src, Size(0, 0), blurRadius, 0, BORDER_DEFAULT );

	LaplacianZC laplacian;
	laplacian.setAperture(kernelRadius);
	laplacian.computeLaplacian(src);
	Mat laplace = laplacian.getLaplacianImage();

	img = mat2Image(laplace, getColorOrderOfRGB32Format(), QImage::Format_RGB888);
#if 0
	  Mat src_gray, dst;
//  int kernel_size = 3;
  int kernel_size = kernelRadius;
  int scale = 1;
  int delta = 0;
  int ddepth = CV_16S;
  
  Mat src = QImageToCvMat(img);

  /// Remove noise by blurring with a Gaussian filter
  //GaussianBlur( src, src, Size(5,5), 0, 0, BORDER_DEFAULT );
  GaussianBlur( src, src, Size(blurRadius, blurRadius), 0, 0, BORDER_DEFAULT );

  /// Convert the image to grayscale
  cvtColor( src, src_gray, COLOR_RGB2GRAY );


  Mat abs_dst;

  Laplacian( src_gray, dst, ddepth, kernel_size, scale, delta, BORDER_DEFAULT );
  convertScaleAbs( dst, abs_dst );

  //img.convertToFormat(QImage::Format_Indexed8);
//img = cvMatToQImage(abs_dst);
  img = mat2Image(abs_dst, getColorOrderOfRGB32Format(), QImage::Format_RGB888);

//img = tmp;
  /// Show what you got
 // imshow( window_name, abs_dst );

 // waitKey(0);
#endif
}

void Filters::zeroCrossing(QImage &img, float blurRadius, float kernelRadius)
{
	Mat src = image2Mat(img, CV_8UC1, getColorOrderOfRGB32Format());
	GaussianBlur( src, src, Size(0, 0), blurRadius, 0, BORDER_DEFAULT );
	//GaussianBlur( src, src, Size(blurRadius, blurRadius), 0, 0, BORDER_DEFAULT );

	LaplacianZC laplacian;
	laplacian.setAperture(kernelRadius);
	Mat flap = laplacian.computeLaplacian(src);

	double mn, mx;
	minMaxLoc(flap, &mn, &mx);

	Mat laplace = laplacian.getLaplacianImage();
	Mat zeros;
	zeros = laplacian.getZeroCrossings(mx); //0	

	img = mat2Image(zeros, getColorOrderOfRGB32Format(), QImage::Format_RGB888);

	/*
		//Mat src = QImageToCvMat(img);;//imread("C:\\qwt-6.1.0\\examples\\imageviewer4\\05-11-2014_13-13-24\\a12-103h_original.jpg", 0);//QImageToCvMat(img).clone();
	Mat src = image2Mat(img, CV_8UC1, getColorOrderOfRGB32Format());
	GaussianBlur( src, src, Size(blurRadius, blurRadius), 0, 0, BORDER_DEFAULT );

	LaplacianZC laplacian;
	laplacian.setAperture(kernelRadius);
	Mat flap = laplacian.computeLaplacian(src);

	double mn, mx;
	minMaxLoc(flap, &mn, &mx);

	Mat laplace = laplacian.getLaplacianImage();
	Mat zeros;
	zeros = laplacian.getZeroCrossings(mx);
	//imshow("lap", laplace);
	img = mat2Image(laplace, getColorOrderOfRGB32Format(), QImage::Format_RGB888);
	*/
#if 0
  Mat src = QImageToCvMat(img), dst, src_gray;
  cvtColor( src, src, COLOR_RGB2GRAY );
  GaussianBlur( src, src_gray, Size(5, 5), 0.5);

  Mat abs_dst;

  Laplacian( src_gray, dst, CV_32F, 3, 1, 0, BORDER_DEFAULT );
  convertScaleAbs( dst, abs_dst );
  imshow("Hell", abs_dst);
  
   float threshold=1.0;
   Mat laplace = dst;

   // Create the iterators
   cv::Mat_<float>::const_iterator it= laplace.begin<float>()+laplace.step1();
   cv::Mat_<float>::const_iterator itend= laplace.end<float>();
   cv::Mat_<float>::const_iterator itup= laplace.begin<float>();

        // Binary image initialize to white
        cv::Mat binary(laplace.size(),CV_8U,cv::Scalar(255));
        cv::Mat_<uchar>::iterator itout= binary.begin<uchar>()+binary.step1();

        // negate the input threshold value
        threshold *= -1.0;

        for ( ; it!= itend; ++it, ++itup, ++itout) {

        // if the product of two adjascent pixel is negative
        // then there is a sign change
        if (*it * *(it-1) < threshold)
        *itout= 0; // horizontal zero-crossing
        else if (*it * *itup < threshold)
        *itout= 0; // vertical zero-crossing
        }
        
        //return ;
		img = mat2Image(binary, getColorOrderOfRGB32Format(), QImage::Format_RGB888);


#if 0
	Mat laplacian = dst;//= QImageToCvMat(img);
	Mat zero_crossings;

	Mat* result = new Mat( laplacian.size(), CV_8U, Scalar(0) );
	zero_crossings = *result;
	int image_rows = laplacian.rows;
	int image_channels = laplacian.channels();
	int values_on_each_row = laplacian.cols * image_channels;
	float laplacian_threshold = 0.0;
	// Find Zero Crossings
	for (int row=1; row < image_rows; row++) {
		float* prev_row_pixel = laplacian.ptr<float>(row-1) +1;
		float* curr_row_pixel = laplacian.ptr<float>(row);
		uchar* output_pixel = zero_crossings.ptr<uchar>(row) +1;
		for (int column=1; column < values_on_each_row; column++)
		{
			float prev_value_on_row = *curr_row_pixel;
			curr_row_pixel++;
			float curr_value = *curr_row_pixel;
			float prev_value_on_column = *prev_row_pixel;
			float difference = 0.0;
			if (((curr_value > 0) && (prev_value_on_row < 0)) ||
				((curr_value < 0) && (prev_value_on_row > 0)))
				difference = abs(curr_value - prev_value_on_row);
			if ((((curr_value > 0) && (prev_value_on_column < 0)) ||
				 ((curr_value < 0) && (prev_value_on_column > 0))) &&
				(abs(curr_value - prev_value_on_column) > difference))
				difference = abs(curr_value - prev_value_on_column);
 			*output_pixel = (difference > laplacian_threshold) ? 255 : 0;// (int) ((100 * difference) / laplacian_threshold);
			prev_row_pixel++;
			output_pixel++;
		}
	}

	img = mat2Image(zero_crossings, getColorOrderOfRGB32Format(), QImage::Format_RGB888);
#endif

#endif
}

void Filters::multiply(QImage &orig, QImage &img1)
{
	int width = img1.width();
	int height = img1.height();


	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			QRgb tmp = orig.pixel(i, j);

			if(QColor(img1.pixel(i, j)).red() == 0)
				img1.setPixel(i, j, 0);
			else
				img1.setPixel(i, j, tmp);


			//orig.setPixel(i, j, tmp);
			/*if(QColor(img1.pixel(i, j)).red() == 0)
				img1.setPixel(i, j, 0);
			else
				img1.setPixel(i, j, tmp);*/
			//int res = QColor(orig.pixel(i, j)).red() - QColor(img1.pixel(i, j)).red();
			//img1.setPixel(i, j, QColor(res, res, res).rgba());
		}
}

void Filters::nullify(QImage &img1)
{
	int pixels[] = {1, 3, 5, 7};
	//int pixels[] = {1, 4, 6, 8, 10};
	//int pixels[] = {5, 7, 3, 1};

	int width = img1.width();
	int height = img1.height();


	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			//QRgb tmp = orig.pixel(i, j);
			for(int z = 0; z < sizeof(pixels)/sizeof(pixels[0]); z++)
			{
				if(QColor(img1.pixel(i, j)).red() == pixels[z])
					img1.setPixel(i, j, 0);
			}
		}
}

void Filters::combine(QImage &orig, QImage &img1)
{
	int width = img1.width();
	int height = img1.height();


	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
		{
			QRgb tmp = orig.pixel(i, j);

			if(QColor(img1.pixel(i, j)).red() == 255)
				img1.setPixel(i, j, tmp);
			
		}
}


#if 0

void Filters::Watershed(QImage &imgTemp, QString &path)
{
	int sz = 2; 
	double C = 0.02;
	//путь к картинке
	QString picturePath = QFileDialog::getOpenFileName(0,
		tr("Open Picture"), QDir::currentPath());

	//путь к маркеру маркера
	QString maskPath = QFileDialog::getOpenFileName(0,
		tr("Open Mask"), QDir::currentPath());
	if (maskPath.isEmpty()) 
		return;

	//загрузка маркера
	IplImage *i1 = cvCreateImage(cvSize(imgTemp.width(), imgTemp.height()), IPL_DEPTH_8U, 1);
	i1 = cvLoadImage(maskPath.toLocal8Bit().constData(), 0);
	IplImage *i2 = cvCreateImage(cvGetSize(i1), IPL_DEPTH_32S, 1);
	cvConvertScale(i1, i2);
	cv::Mat mat2(i2);
	CvMat mmm2 = mat2;

	//повторная загрузка картинки
	IplImage *i0 = cvLoadImage(picturePath.toLocal8Bit().constData(), 1);
	cv::Mat mat(i0);

	CvMat mmm = mat;


	IplImage *clone = cvLoadImage(picturePath.toLocal8Bit().constData(), 1);
	cvZero(clone);

	const int IN_QUEUE = -2;
	const int WSHED = -1;
	const int NQ = 256;

	const int thres = 256;
	const int block = (2*sz+1)*(2*sz+1);
	CvMemStorage* storage = 0;

	CvMat sstub, *src;
	CvMat dstub, *dst;
	CvSize size;
	CvWSNode* free_node = 0, *node;
	CvWSQueue q[NQ];
	int active_queue;
	int i, j;
	int db, dg, dr;
	int* mask;
	uchar* img;
	int mstep, istep;
	int subs_tab[513];


#define ws_max(a,b) ((b) + subs_tab[(a)-(b)+NQ])
#define ws_min(a,b) ((a) - subs_tab[(a)-(b)+NQ])
#define ws_push(idx,mofs,iofs)  \
	{                               \
	if( !free_node )            \
	free_node = icvAllocWSNodes( storage );\
	node = free_node;           \
	free_node = free_node->next;\
	node->next = 0;             \
	node->mask_ofs = mofs;      \
	node->img_ofs = iofs;       \
	if( q[idx].last )           \
	q[idx].last->next=node; \
	else                        \
	q[idx].first = node;    \
	q[idx].last = node;         \
	}

#define ws_pop(idx,mofs,iofs)   \
	{                               \
	node = q[idx].first;        \
	q[idx].first = node->next;  \
	if( !node->next )           \
	q[idx].last = 0;        \
	node->next = free_node;     \
	free_node = node;           \
	mofs = node->mask_ofs;      \
	iofs = node->img_ofs;       \
	}

#define c_diff(ptr1,ptr2,diff)      \
	{                                   \
	db = abs((ptr1)[0] - (ptr2)[0]);\
	dg = abs((ptr1)[1] - (ptr2)[1]);\
	dr = abs((ptr1)[2] - (ptr2)[2]);\
	diff = ws_max(db,dg);           \
	diff = ws_max(diff,dr);         \
	assert( 0 <= diff && diff <= 255 ); \
	}


	src = &mmm;
	dst = &mmm2;

	size = cvGetMatSize(src);

	storage = cvCreateMemStorage();

	istep = src->step;
	img = src->data.ptr;
	mstep = dst->step / sizeof(mask[0]);
	mask = dst->data.i;

	memset( q, 0, NQ*sizeof(q[0]) );

	for( i = 0; i < 256; i++ )
		subs_tab[i] = 0;
	for( i = 256; i <= 512; i++ )
		subs_tab[i] = i - 256;

	//рисование рамки из пикселей вокруг изображения (border of dummy "watershed")
	//записывается число -1 (255). WSHED - константа.
	for( i = 0; i < sz; i++ ) 
		for( j = 0; j < size.width; j++ ){
			mask[j+i*mstep] = mask[j + mstep*(size.height-i-1)] = WSHED;
		}

		mask += (sz-1)*mstep;
		img += (sz-1)*istep;

		//начальная фаза. все соседние пиксели каждого маркера помещаются в упорядоченную очередь
		//определяются начальные границы бассейнов
		for( i = sz; i < size.height-sz; i++ )
		{
			img += istep; mask += mstep;
			for( int n=0; n < sz; n++ )
				mask[0+n] = mask[size.width-n-1] = WSHED;

			for( j = sz; j < size.width-sz; j++ )
			{
				int* m = mask + j;
				if( m[0] < 0 ) m[0] = 0;
				if( m[0] == 0 && (m[-1] > 0 || m[1] > 0 || m[-mstep] > 0 || m[mstep] > 0) )
				{
					uchar* ptr = img + j*3;
					int idx = 256, t;
					if( m[-1] > 0 )
						c_diff( ptr, ptr - 3, idx );
					if( m[1] > 0 )
					{
						c_diff( ptr, ptr + 3, t );
						idx = ws_min( idx, t );
					}
					if( m[-mstep] > 0 )
					{
						c_diff( ptr, ptr - istep, t );
						idx = ws_min( idx, t );
					}
					if( m[mstep] > 0 )
					{
						c_diff( ptr, ptr + istep, t );
						idx = ws_min( idx, t );
					}
					assert( 0 <= idx && idx <= 255 );
					ws_push( idx, i*mstep + j, i*istep + j*3 );
					m[0] = IN_QUEUE;
				}
			}
		}

		//поиск первой непустой очереди
		for( i = 0; i < NQ; i++ )
			if( q[i].first )
				break;

		//если нет маркеров, то выход, т.е. :
		//число очередей - 255, если NQ будет 256, то выход
		if( i == NQ )
			exit;

		active_queue = i;
		img = src->data.ptr;
		mask = dst->data.i;

		//рекурсивное заполнение бассейнов
		for(;;)
		{
			int mofs, iofs;
			int lab = 0, t, t1, b;
			int* m;
			uchar* ptr;

			if( q[active_queue].first == 0 )
			{
				for( i = active_queue+1; i < NQ; i++ )
					if( q[i].first )
						break;

				if( i == NQ )
					break;
				active_queue = i;
			}

			ws_pop( active_queue, mofs, iofs );

			m = mask + mofs;
			ptr = img + iofs;
			t = m[-1];
			if( t > 0 ) lab = t;
			t = m[1];
			if( t > 0 )
				if( lab == 0 ) lab = t;
				else if( t != lab ) lab = WSHED;
				t = m[-mstep];
				if( t > 0 )
					if( lab == 0 ) lab = t;
					else if( t != lab ) lab = WSHED;
					t = m[mstep];
					if( t > 0 )
						if( lab == 0 ) lab = t;
						else if( t != lab ) lab = WSHED;

						assert( lab != 0 );
						m[0] = lab;
						if( lab == WSHED )
							continue;

						if( m[-1] == 0 )
						{
							c_diff( ptr, ptr - 3, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr-3, ptr-i*3+j*istep, t1);
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);   
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs - 1, iofs - 3 );
							active_queue = ws_min( active_queue, t );
							m[-1] = IN_QUEUE;
						}
						if( m[1] == 0 )
						{
							c_diff( ptr, ptr + 3, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr+3, ptr+i*3+j*istep, t1);
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs + 1, iofs + 3 );
							active_queue = ws_min( active_queue, t );
							m[1] = IN_QUEUE;
						}
						if( m[-mstep] == 0 )
						{
							c_diff( ptr, ptr-istep, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr-istep, ptr+j*3-i*istep, t1 );
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs - mstep, iofs - istep );
							active_queue = ws_min( active_queue, t );
							m[-mstep] = IN_QUEUE;
						}
						if( m[mstep] == 0 )
						{
							c_diff( ptr, ptr+istep, t );
							if( t < thres ){
								double a = 0;
								for( int i=-sz; i<=sz; i++ )
									for( int j=-sz; j<=sz; j++ ){
										c_diff( ptr+istep, ptr+j*3+i*istep, t1 );
										a = a + exp(-C*double(t1));
									}
									b = int(log(double(block)-a+1.0)/C+0.5);
									t=t+b;
							}
							if( t>255 ) t=255;
							ws_push( t, mofs + mstep, iofs + istep );
							active_queue = ws_min( active_queue, t );
							m[mstep] = IN_QUEUE;
						}
		}
		cvReleaseMemStorage( &storage );

		//получение результирующего изображения 
		//обработанного водоразделом
		for( int i = 0; i < i2->height; i++ )
			for( int j = 0; j < i2->width; j++ )
			{
				int idx = CV_IMAGE_ELEM( i2, int, i, j ); //i2 - финальное изображение с бассейнами водоразделов
				uchar* dst1 = &CV_IMAGE_ELEM( clone, uchar, i, j*3 );

				if( idx == -1 )
				{
					dst1[0] = (uchar)0; dst1[1] = (uchar)255; dst1 [2] = (uchar)0;
					if (i+1< i2->height && i-1>=0 && j+1<i2->width && j-1>=0) 
					{
						for (int n=-1; n<=1; n++)
							for (int m=-1; m<=1; m++)
							{
								uchar* dst2 = &CV_IMAGE_ELEM( clone, uchar, i+n, (j+m)*3 );
								dst2[0] = (uchar)0; dst2[1] = (uchar)0; dst2[2] = (uchar)255;
								//отмечаем красную рамочку на исходном изображении
								imgTemp.setPixel(j+m, i+n, QColor(255, 0, 0).rgba());
							}
					}
				}
			}
			/*
		for( int i = 10; i < i2->height-10; i++ )
			for( int j = 10; j < i2->width-10;  )
			{
				if(QColor(255, 0, 0) == imgTemp.pixel(j, i))
				{
					for( ; j < i2->width-10; j++ )
					{
						imgTemp.setPixel(j, i, QColor(0, 0, 0).rgba());
						if(QColor(255, 0, 0) == imgTemp.pixel(j, i))
							break;
					}
				}
				else
					j++;
			}
		*/
		cvConvertScale(i2, i1, 1, 0);
		cvSaveImage("final_markers.png", i1);

}
#endif