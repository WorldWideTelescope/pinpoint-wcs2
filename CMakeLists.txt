cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

set(VERSION "1.0.0")
project(PinpointWCS VERSION ${VERSION})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt5 COMPONENTS Core Network Widgets REQUIRED)

find_package(CFITSIO REQUIRED)

find_package(WCSLIB REQUIRED)

# Yes, libwcs and wcslib are different things, and we currently need both of
# them.
find_package(LIBWCS REQUIRED)

find_package(Eigen3 REQUIRED)

find_package(XPA REQUIRED)

include_directories(
    ${CMAKE_SOURCE_DIR}/backend
    ${CMAKE_SOURCE_DIR}/gui
    ${CMAKE_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${WCSLIB_INCLUDES}
    ${LIBWCS_INCLUDE_DIR}
    ${CFITSIO_INCLUDE_DIR}
    ${XPA_INCLUDE_DIR}
)

add_executable(
    PinpointWCS
    PinpointWCS.qrc
    main.cpp
    backend/ComputeWCS.cpp
    backend/CoordinateDelegate.cpp
    backend/CoordinateModel.cpp
    backend/DS9Thread.cpp
    backend/EpoImage.cpp
    backend/ExportWCS.cpp
    backend/FitsImage.cpp
    backend/FITSThread.cpp
    backend/PinpointWCSUtils.cpp
    backend/PPWcsImage.cpp
    backend/RemoteData.cpp
    backend/ComputeWCS.h
    backend/CoordinateDelegate.h
    backend/CoordinateModel.h
    backend/DS9Thread.h
    backend/EpoImage.h
    backend/ExportWCS.h
    backend/FitsImage.h
    backend/FITSThread.h
    backend/PinpointWCSUtils.h
    backend/PPWcsImage.h
    backend/RemoteData.h
    gui/AboutDialog.ui
    gui/CoordinatePanel.ui
    gui/CoordinateTableDialog.ui
    gui/FitsToolbar.ui
    gui/HelpPanel.ui
    gui/PinpointWCS.ui
    gui/WcsInfoPanel.ui
    gui/AboutDialog.cpp
    gui/Commands.cpp
    gui/CoordinateMarker.cpp
    gui/CoordinatePanel.cpp
    gui/CoordinateTableDialog.cpp
    gui/DropArea.cpp
    gui/FitsToolbar.cpp
    gui/GraphicsScene.cpp
    gui/GraphicsView.cpp
    gui/HelpPanel.cpp
    gui/mainwindow.cpp
    gui/MessageBox.cpp
    gui/WcsInfoPanel.cpp
    gui/AboutDialog.h
    gui/Commands.h
    gui/CoordinateMarker.h
    gui/CoordinatePanel.h
    gui/CoordinateTableDialog.h
    gui/DropArea.h
    gui/FitsToolbar.h
    gui/GraphicsScene.h
    gui/GraphicsView.h
    gui/HelpPanel.h
    gui/mainwindow.h
    gui/MessageBox.h
    gui/WcsInfoPanel.h
)

target_link_libraries(
    PinpointWCS
    PUBLIC
    Qt5::Network
    Qt5::Widgets
)
