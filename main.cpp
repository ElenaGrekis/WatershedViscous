
#include <QApplication>

#include "mainwindow.h"

void myMessageOutput(QtMsgType type, const char *msg)
{}

int main(int argc, char *argv[])
{
	qInstallMsgHandler(myMessageOutput);
	QApplication app(argc, argv);

	QTextCodec::setCodecForCStrings ( QTextCodec::codecForName("Windows-1251") );
	QTextCodec::setCodecForLocale ( QTextCodec::codecForName("Windows-1251") );
    

	MainWindow m;
	m.show();
    return app.exec();
}
