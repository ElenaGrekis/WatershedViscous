#include "canvas.h"


Plot::Plot( QWidget *parent, const QString &  name):
    QwtPlot( parent )
{
    setAutoFillBackground( true );
    //setPalette( QPalette( QColor( 165, 193, 228 ) ) );
    //updateGradient();

    setTitle( name );
    insertLegend( new QwtLegend(), QwtPlot::RightLegend );

    // axes
    setAxisTitle( xBottom, "x -->" );
    setAxisScale( xBottom, 0.0, 10.0 );

    setAxisTitle( yLeft, "y -->" );
    setAxisScale( yLeft, -1.0, 1.0 );

    // canvas
    QwtPlotCanvas *canvas = new QwtPlotCanvas();
    canvas->setLineWidth( 1 );
    canvas->setFrameStyle( QFrame::Box | QFrame::Plain );
    canvas->setBorderRadius( 15 );

    //QPalette canvasPalette( Qt::white );
    //canvasPalette.setColor( QPalette::Foreground, QColor( 133, 190, 232 ) );
    //canvas->setPalette( canvasPalette );

    setCanvas( canvas );

    // panning with the left mouse button
    ( void ) new QwtPlotPanner( canvas );

    // zoom in/out with the wheel
    ( void ) new QwtPlotMagnifier( canvas );

  //  populate();
}

void Plot::updateGradient()
{
    QPalette pal = palette();

    const QColor buttonColor = pal.color( QPalette::Button );

    QLinearGradient gradient( rect().topLeft(), rect().bottomLeft() );
    gradient.setColorAt( 0.0, Qt::white );
    gradient.setColorAt( 0.7, buttonColor );
    gradient.setColorAt( 1.0, buttonColor );

    pal.setBrush( QPalette::Window, gradient );
    setPalette( pal );
}