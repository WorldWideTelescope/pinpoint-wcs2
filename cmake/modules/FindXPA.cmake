# Find the XPA library
#
# This module defines these variables:
#
#   XPA_FOUND
#      True if the XPA library was found.
#   XPA_LIBRARY
#      The location of the XPA library.
#   XPA_INCLUDE_DIR
#      The include path of the XPA library.

#
# Find the header file
#
FIND_PATH(XPA_INCLUDE_DIR xpa.h)

#
# Find the library
#
FIND_LIBRARY(XPA_LIBRARY xpa)

SET(XPA_FOUND false)
IF(XPA_INCLUDE_DIR AND XPA_LIBRARY)
   SET(XPA_FOUND true)
ENDIF(XPA_INCLUDE_DIR AND XPA_LIBRARY)

IF(XPA_FOUND)
   MESSAGE(STATUS "Found XPA: ${XPA_LIBRARY}")
ELSE(XPA_FOUND)
   MESSAGE(FATAL_ERROR "Could not find the XPA library")
ENDIF(XPA_FOUND)
