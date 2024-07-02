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

#ifndef COORDINATE_MARKER_H
#define COORDINATE_MARKER_H

#include <QGraphicsItem>
#include <QStyleOptionGraphicsItem>

QT_BEGIN_NAMESPACE
class QPersistentModelIndex;
QT_END_NAMESPACE

class CoordinateMarker : public QGraphicsItem
{
	
	friend class MoveCommand;
	friend class GraphicsScene;
	
public:
	enum {Type = UserType + 1};
    CoordinateMarker(QModelIndex &idx, QGraphicsItem *parent = nullptr);
    ~CoordinateMarker() override;
	
	// Required methods to implement
    QRectF boundingRect() const override;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *item, QWidget *widget) override;
    QPainterPath shape() const override;
    int type() const override { return Type; } //enables qgraphicsitem_cast

protected:
    void wheelEvent(QGraphicsSceneWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override;
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override;
	
private:
	QPersistentModelIndex *index;
	void setRadius();
	float setPenWidth();
	float radius;
	float penWidth;
	float measure;
	float scale;
	
};

#endif
