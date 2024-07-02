# PinpointWCS

## Building

We now use CMake. My personal CMake command:

```
cmake \
  -DCFITSIO_INCLUDE_DIR=/a/dasch/pfx/include -DCFITSIO_LIBRARY=/a/dasch/pfx/lib/libcfitsio.a \
  -DLIBWCS_INCLUDE_DIR=/a/dasch/pfx/include -DLIBWCS_LIBRARY=/a/dasch/pfx/lib/libwcs.a \
  -DEIGEN3_INCLUDE_DIR=/a/include/eigen3 \
  -DXPA_INCLUDE_DIR=/a/include -DXPA_LIBRARY=/a/lib/libxpa.a
```

Dependencies:

- [Qt 5](https://doc.qt.io/qt-5/)
- [expat](https://libexpat.github.io/)
- [libwcs](http://tdc-www.harvard.edu/wcstools/subroutines/libwcs.wcs.html) from [wcstools](http://tdc-www.harvard.edu/wcstools/)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/)
- [Eigen 3](https://eigen.tuxfamily.org/)
- [XPA](https://github.com/ericmandel/xpa)
- [zlib](http://www.zlib.net/)
- We vendor in the [Adobe XMP Toolkit
  SDK](https://github.com/adobe/XMP-Toolkit-SDK/)