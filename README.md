# ghost_fluid_method
This repository contain the codes (written in C++) for the Ghost Fluid Method (GFM) of Liu et al. (1999, JCP). Finite Difference Method is used to discretize the Poisson equation.

The codes (ps_1a.cpp, ps_1b.cpp, ps_2.cpp, etc.) are designed to solve the test problems (Examples 1, 2, ..., resectively) of Liu et al. (1999, JCP).

# Reference
Liu, X.D., Fedkiw, R.P. and Kang, M., 2000. A boundary condition capturing method for Poisson's equation on irregular domains. Journal of computational Physics, 160(1), pp.151-178.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

The process is similar for all the codes (ps_*.cpp). To compile the code ps_1a.cpp

```bash
g++ ps_1a.cpp -o output
```
To run this code

```bash
./output
```
