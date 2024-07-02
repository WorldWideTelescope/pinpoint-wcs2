# PinpointWCS

## Building

We now use CMake. My personal CMake command:

```
cmake \
  -DCFITSIO_INCLUDE_DIR=/a/include -DCFITSIO_LIBRARY=/a/lib/libcfitsio.a \
  -DWCSLIB_ROOT_DIR=/a \
  -DEIGEN3_INCLUDE_DIR=/a/include/eigen3 \
  -DXPA_INCLUDE_DIR=/a/include -DXPA_LIBRARY=/a/lib/libxpa.a
```
