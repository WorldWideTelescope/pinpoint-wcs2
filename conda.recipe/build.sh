#! /bin/bash

set -xeuo pipefail

# Workaround/hack -- the pkgw-forge version of wcstools puts its headers in
# $PREFIX/include/wcstools/libwcs/*.h, while the conda-forge package does
# not include the "libwcs/" component. The upstream source currently uses
# #include statements of the form `#include <libwcs/wcs.h>`. So just copy
# files around to make everything appear homogeneous, if needed.
if [ ! -d "$PREFIX/include/wcstools/libwcs" ] ; then
    mkdir "$PREFIX/include/wcstools/libwcs"
    cp -av "$PREFIX/include/wcstools"/*.h "$PREFIX/include/wcstools/libwcs/"
fi

cmake_args=(
    -DCFITSIO_INCLUDE_DIR=$PREFIX/include
    -DCFITSIO_LIBRARY=$PREFIX/lib/libcfitsio$SHLIB_EXT
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_COLOR_MAKEFILE=OFF
    -DCMAKE_INSTALL_PREFIX=$PREFIX
    -DEIGEN3_INCLUDE_DIR=$PREFIX/include/eigen3
    -DLIBWCS_INCLUDE_DIR=$PREFIX/include/wcstools
    -DLIBWCS_LIBRARY=$PREFIX/lib/libwcstools$SHLIB_EXT
    -DXPA_INCLUDE_DIR=$PREFIX/include
    -DXPA_LIBRARY=$PREFIX/lib/libxpa.a
)

if [[ $(uname) == Darwin ]] ; then
    linkflags="-Wl,-rpath,$PREFIX/lib $LDFLAGS"

    # cmake_args+=(
    #     -DCMAKE_CXX_FLAGS="-arch $OSX_ARCH -stdlib=libc++ -D_LIBCPP_DISABLE_AVAILABILITY"
    #     -DCMAKE_OSX_DEPLOYMENT_TARGET=$MACOSX_DEPLOYMENT_TARGET
    #     -DCMAKE_OSX_SYSROOT=/
    # )
else
    linkflags="-Wl,-rpath-link,$PREFIX/lib -ldl $LDFLAGS"
fi

cmake_args+=(
    -DCMAKE_CXX_FLAGS="$CXXFLAGS -I$(pwd)/build/adobexmp"
    -DCMAKE_EXE_LINKER_FLAGS="$linkflags"
    -DCMAKE_MODULE_LINKER_FLAGS="$linkflags"
    -DCMAKE_SHARED_LINKER_FLAGS="$linkflags"
)

# temp hack for my dumb tarball
rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt .git

mkdir build
cd build

cmake "${cmake_args[@]}" ..
make -j${CPU_COUNT} VERBOSE=1

mkdir -p $PREFIX/bin
cp -p PinpointWCS $PREFIX/bin

# Avoid packaging these files if we created them above
rm -rf "$PREFIX/include/wcstools/libwcs"
