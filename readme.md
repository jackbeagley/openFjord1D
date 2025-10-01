# openFjord1D
This is a program that models the evolution of a fjord basin.
The 1-dimensional equation for diffusion is used with a term that deals with the changing horizontal cross-sectional area of a basin.
Salinity and temperature are specified at the top of the silled-basin and evolve inside the basin, determined by a specified diffusivity.
Density is calculated by https://github.com/TEOS-10/GSW-Fortran.
If water at sill depth is denser than water inside of the basin, lighter water in the basin is replace by water of the same properties of the sill water.

## Building
Use the following commands to build.
A Fortran compiler is required, as well as python
```
cmake -S src -B build
cmake --build build
```

If everything goes to plan, this should generate a compiled library within the build directory (e.g. `build/openFjord1D.cpython-313-darwin.so` on MacOS).
This library can be copied to the location of the script that is using it and then imported by that script:

```python
import openFjord1D as of
```
