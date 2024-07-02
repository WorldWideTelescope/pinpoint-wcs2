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
