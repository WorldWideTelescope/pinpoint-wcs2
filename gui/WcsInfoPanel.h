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

#ifndef WCS_INFO_PANEL_H
#define WCS_INFO_PANEL_H

#include <QFrame>
#include "ui_WcsInfoPanel.h"
#include "wcs.h"
#include "fitsfile.h"


class WcsInfoPanel : public QFrame
	{
		Q_OBJECT
		
	public:
        WcsInfoPanel(bool ref, QWidget *parent = nullptr);
        ~WcsInfoPanel() override;
        void loadWCS(struct WorldCoor* wcs, double rms_x = 0.0, double rms_y = NULL);
		void clear();
		
	public slots:
		void parentResized(QSize sz);
		
	private:
		// Attributes
		Ui::WcsInfoPanel ui;
		bool reference;
		
		// Methods
		void updateFontSize(QFont font);
	};

#endif
