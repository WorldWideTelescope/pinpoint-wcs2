/*
 *  PinpointWCS is developed by the Chandra X-ray Center
 *  Education and Public Outreach Group
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef GRAPHICSVIEW_H
#define GRAPHICSVIEW_H

#include <QGraphicsView>
#include <QGraphicsItem>
#include <QKeyEvent>

QT_BEGIN_NAMESPACE
class QPixmap;
class QGraphicsPixmapItem;
QT_END_NAMESPACE

#define MINZOOM 0.07
#define MAXZOOM 400

class GraphicsView : public QGraphicsView
{
	Q_OBJECT
	
public:
    GraphicsView(QWidget *parent = nullptr);
    ~GraphicsView() override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;
    virtual void enterEvent(QKeyEvent *event);
    void leaveEvent(QEvent *event) override;
    qreal scaling();
	
	// Attributes
	int rotateFactor;
	
public slots:
	void rotateCW();
	void rotateCCW();

protected:
	// Methods
    void mouseReleaseEvent(QMouseEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
	void scaleView(qreal scaleFactor);
	
signals:
	void objectResized(QSize s);
	void mouseEnterEvent(GraphicsView *gv);

};

#endif
