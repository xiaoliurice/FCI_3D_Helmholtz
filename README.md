# 3D_Helmholtz_solver

Fast contour integration (FCI) method for solving the Helmholtz equation

Developers: Xiao Liu, Yuanzhe Xi

Consider the linear system (A-M)u=f, where A is the discrete Laplacian (here using FFT or finite difference), M is the potential (kh^2). The contour integration approach combines the solutions of multiple shifted systems (A-zM)u = f. The imaginary part of z makes the shifted systems easy to solve. A fixed-point iteration is used to solve shifted problems.

This repository contains Matlab and C code.
