cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(VERSION "1.0.0")
project(PinpointWCS VERSION ${VERSION})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
add_definitions(${Qt5Core_DEFINITIONS})
include_directories(${Qt5Core_INCLUDE_DIRS})

include_directories(
	${CMAKE_SOURCE_DIR}/gui
)

add_executable(
    PinpointWCS
    PinpointWCS.qrc
    main.cpp
    gui/AboutDialog.ui
    gui/CoordinatePanel.ui
    gui/CoordinateTableDialog.ui
    gui/FitsToolbar.ui
    gui/HelpPanel.ui
    gui/PinpointWCS.ui
    gui/WcsInfoPanel.ui
    gui/mainwindow.h
)

target_link_libraries(
    PinpointWCS
	PUBLIC
	Qt5::Widgets
)
