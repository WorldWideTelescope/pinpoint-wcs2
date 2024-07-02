CONFIG += release
#CONFIG += debug_and_release

macx {

        INCLUDEPATH += . gui backend ../cfitsio/include ../libwcs ../eigen
	INCLUDEPATH += ../XMP-Toolkit-SDK-main/public/include
	INCLUDEPATH += ../xpa-master/include

	LIBS += -L../cfitsio/lib -lcfitsio
	LIBS += -L../libwcs -lwcs -lm
	LIBS += -L../XMP-Toolkit-SDK-main/public/libraries/macintosh/intel_64_libcpp/Debug
	LIBS += -lXMPCoreStatic -lXMPFilesStatic
	LIBS += -framework CoreServices
	LIBS += -L../xpa-master/lib -lxpa

	DEFINES += MAC_ENV="1"
}

unix:!macx {

	INCLUDEPATH += . gui backend ../cfitsio/include ../libwcs/ ../eigen
	INCLUDEPATH += ../XMP-Toolkit-SDK-main/public/include
	INCLUDEPATH += ../xpa-master/include

	LIBS += -L../cfitsio/lib -lcfitsio
	LIBS += -L../libwcs -lwcs -lm
	LIBS += ../XMP-Toolkit-SDK-main/public/libraries/intel_64_libcpp/Debug/staticXMPCore.ar
	LIBS += ../XMP-Toolkit-SDK-main/public/libraries/intel_64_libcpp/Debug/staticXMPFiles.ar
	LIBS += -L../xpa-master/lib -lxpa

	DEFINES += UNIX_ENV="1"
}

win32 {

	INCLUDEPATH += . gui backend ../cfitsiodll_3280_vcc ../libwcs/ ../eigen
	INCLUDEPATH += ../XMP-Toolkit-SDK-main/public/include
	INCLUDEPATH += ../xpa-master/include

	LIBS += -L../cfitsiodll_3280_vcc -lcfitsio
	LIBS += -L../libwcs -lwcs -lm
	LIBS += -L../XMP-Toolkit-SDK-5.1.2/public/libraries/windows/release
	LIBS += -lXMPCoreStaticRelease -lXMPFilesStaticRelease
	LIBS += -L../xpa-2.1.13/lib -lxpa
	
	DEFINES += WIN_ENV="1"
}