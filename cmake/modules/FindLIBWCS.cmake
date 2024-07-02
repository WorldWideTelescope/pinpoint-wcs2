# Find the LIBWCS library
#
# This module defines these variables:
#
#   LIBWCS_FOUND
#      True if the LIBWCS library was found.
#   LIBWCS_LIBRARY
#      The location of the LIBWCS library.
#   LIBWCS_INCLUDE_DIR
#      The include path of the LIBWCS library.

#
# Find the header file
#
FIND_PATH(LIBWCS_INCLUDE_DIR libwcs/fitsfile.h)

#
# Find the library
#
FIND_LIBRARY(LIBWCS_LIBRARY wcs)

SET(LIBWCS_FOUND false)
IF(LIBWCS_INCLUDE_DIR AND LIBWCS_LIBRARY)
   SET(LIBWCS_FOUND true)
ENDIF(LIBWCS_INCLUDE_DIR AND LIBWCS_LIBRARY)

IF(LIBWCS_FOUND)
   MESSAGE(STATUS "Found LIBWCS: ${LIBWCS_LIBRARY}")
ELSE(LIBWCS_FOUND)
   MESSAGE(FATAL_ERROR "Could not find the LIBWCS library")
ENDIF(LIBWCS_FOUND)
