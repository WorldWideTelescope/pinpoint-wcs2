######################################################################
# Automatically generated by qmake (2.01a) Mon Oct 25 16:26:12 2010
######################################################################

TEMPLATE = app
TARGET =
DEPENDPATH += . backend gui
INCLUDEPATH += . gui backend ../cfitsio/include ../libwcs/ ../eigen
INCLUDEPATH += ../XMP-Toolkit-SDK-5.1.2/public/include
CONFIG += release
#CONFIG += debug_and_release
LIBS += -L../cfitsio/lib -lcfitsio
LIBS += -L../libwcs -lwcs -lm
LIBS += -L../XMP-Toolkit-SDK-5.1.2/public/libraries/macintosh/release
LIBS += -lXMPCoreStaticRelease -lXMPFilesStaticRelease
LIBS += -framework CoreServices
DEFINES +=  MAC_ENV="1"
QT += network

# Input
HEADERS += version.h \
           backend/ComputeWCS.h \
           backend/CoordinateDelegate.h \
           backend/CoordinateModel.h \
           backend/EpoImage.h \
           backend/ExportWCS.h \
           backend/FitsImage.h \
           backend/PinpointWCSUtils.h \
           backend/PPWcsImage.h \
           backend/RemoteData.h \
           gui/AboutDialog.h \
           gui/Commands.h \
           gui/CoordinateMarker.h \
           gui/CoordinatePanel.h \
           gui/CoordinateTableDialog.h \
           gui/DropArea.h \
           gui/FitsToolbar.h \
           gui/GraphicsScene.h \
           gui/GraphicsView.h \
           gui/mainwindow.h \
           gui/WcsInfoPanel.h
FORMS += gui/AboutDialog.ui \
         gui/CoordinatePanel.ui \
         gui/CoordinateTableDialog.ui \
         gui/FitsToolbar.ui \
         gui/PinpointWCS.ui \
         gui/WcsInfoPanel.ui
SOURCES += main.cpp \
           backend/ComputeWCS.cpp \
           backend/CoordinateDelegate.cpp \
           backend/CoordinateModel.cpp \
           backend/EpoImage.cpp \
           backend/ExportWCS.cpp \
           backend/FitsImage.cpp \
           backend/PinpointWCSUtils.cpp \
           backend/PPWcsImage.cpp \
           backend/RemoteData.cpp \
           gui/AboutDialog.cpp \
           gui/Commands.cpp \
           gui/CoordinateMarker.cpp \
           gui/CoordinatePanel.cpp \
           gui/CoordinateTableDialog.cpp \
           gui/DropArea.cpp \
           gui/FitsToolbar.cpp \
           gui/GraphicsScene.cpp \
           gui/GraphicsView.cpp \
           gui/mainwindow.cpp \
           gui/WcsInfoPanel.cpp
RESOURCES += PinpointWCS.qrc
