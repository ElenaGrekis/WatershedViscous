#include <QtGui>

#include "appwindow.h"

//! [0]
AppWindow::AppWindow()
{
#if 0
	QWidget *main = new QWidget;
	QWidget *main2 = new QWidget;
	

	vbox = new QVBoxLayout(this);

	hist1 = new Histogram(this);
	hist2 = new Histogram(this);

		
	//vbox->addLayout(hbox2);
	//vbox->addLayout(hbox1);

	main->setLayout(hbox1);
	//main2->setLayout(hbox2);
	//main2->resize(640, 480);
	//main2->show();



	Q_ASSERT(
		connect(view1, 		SIGNAL(sendOriginalImage(QLabel *)), 		
		view2, 		SLOT(recieveImage(QLabel *)) )
		);

	//Q_ASSERT(
		connect(view1, 		SIGNAL(signalSetArray1(int *, int, int)), 		
		hist1, 		SLOT(setArray(int *, int, int)) );
	//	);

	//Q_ASSERT(
		connect(view2, 		SIGNAL(signalSetArray2(int *, int, int)), 		
		hist2, 		SLOT(setArray(int *, int, int)) );
	//	);

    setCentralWidget(main);


    setWindowTitle(tr("Image Viewer"));
    resize(500, 400);
#endif
}
