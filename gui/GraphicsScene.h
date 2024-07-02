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

#ifndef GRAPHICSSCENE_H
#define GRAPHICSSCENE_H

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include "CoordinateMarker.h"
#include "GraphicsView.h"

class GraphicsScene : public QGraphicsScene
{
	Q_OBJECT
		
public:
    GraphicsScene(QPixmap *pix, bool ref, QObject *parent = nullptr);
    ~GraphicsScene() override;
	bool reference;
	float measure;
	QGraphicsRectItem *centralItem;
	bool clickable;
	void signalItemMoved(CoordinateMarker *m, QPointF oldPos);
		
protected:
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event) override;
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

private:
	QPointF oldPos;
	QGraphicsItem *movingItem;
	QGraphicsPixmapItem *ptr_pixmap;
	
	float computeRadii();
	
signals:
	void mousePositionChanged(QPointF pos);
	void sceneDoubleClicked(GraphicsScene *scene, QPointF pos);
	void toggleNeighborScene(bool sendSignal = false);
//    void itemMoved(CoordinateMarker *movedItem, const QPointF &movedFromPosition);
	void itemMoved(GraphicsScene *scene, QPointF newPosition, QPointF oldPosition);
	void markerChange(CoordinateMarker *m);
	void currentSelection(int row);
	void clearCorrespondingSelection();
	void itemPos(QPointF pos);

public slots:
	void updatePixmap(QPixmap *pm);
	void toggleClickable(bool sendSignal = true);
	void findSelectedItem();
	void matchSelectedItem(int row);
	void selectedItemPos();
	
};

#endif
