#pragma once

#include <qapplication.h>
#include <qlayout.h>
#include <qwt_plot.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_point_data.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_magnifier.h>
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <qwt_math.h>
#include <math.h>

class Plot : public QwtPlot
{
public:
    Plot( QWidget *parent = NULL, const QString &  name = "" );

protected:
    //virtual void resizeEvent( QResizeEvent * );

private:
 //   void populate();
    void updateGradient();
};
