# fcimatfree

Fast contour integration preconditioner
Xiao Liu, and Yuanzhe Xi

This version is matrix-free designed for the Helmholtz equation

We consider the linear system (A-M)u=f, where A is discrete Laplacian (here using DST-I), M is the potential (kh^2). The contour integration approach conbines the solutions of multiple shifted systems (A-zM)u = f. The imaginary part of z makes the shifted systems easy to solve. A constant coefficient solver is used as a preconditioner.

