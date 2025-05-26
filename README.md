# DEigen.jl
[![Build Status](https://github.com/ykkan/DEigen.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ykkan/DEigen.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package computes the derivatives of a parametric eigenvalue problem with respect to a single scalar variable.

## Description
Consider a parametric square matrix $A(x)$ that depends on a single scalar variable $x$ and its associated eigenvalue problem:
```math
A(x)v(x) = \lambda(x)v(x). 
```
At a given point $x_0$, if all the derivatives of $A$ to the $n$-th order are known (i.e., $[A^{(k)}(x_0) \mid k=1\ldots K]$), 
then the derivatives of the eigenvalue $[\lambda^{(k)}(x_0) \mid k=1\ldots K]$ and the eigenvector $[v^{(k)}(x_0) \mid k=1\ldots K]$ up to $K$ can be computed. 

## Installation
The package can be installed using Julia's REPL
```julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/ykkan//DEigen.jl.git")
```
or with Pkg mode (hitting `]` in the command prompt)
```julia
pkg> add https://github.com/ykkan//DEigen.jl.git
```

## Usage
### Derivatives of all eigenvalues and eigenvectors
Assume that we want to evaluate the derivatives of all eigenvalues and all eigenvalues for the matrix below at $x=1$ to the 2nd order:
```math
A(x) = 
\begin{bmatrix}
x & 1 \\
1 & 2
\end{bmatrix}
```

```julia
using DEigen
# A0 = A(1),  A1 = A'(1), A2 = A''(1)
A0 = [1.0 1.0;
      1.0 2.0]

A1 = [1.0 0.0;
      0.0 0.0]

A2 = [0.0 0.0;
      0.0 0.0]

values_list, vectors_list = deign(A0, A1, A2)

#=
DEigen{Float64, 2}
values_list:
2×3 Matrix{Float64}:
 0.381966  0.723607  -0.178885
 2.61803   0.276393   0.178885

vectors_list:
2×2×3 Array{Float64, 3}:
[:, :, 1] =
 -0.850651  0.525731
  0.525731  0.850651

[:, :, 2] =
 0.105146   0.17013
 0.17013   -0.105146

[:, :, 3] =
 -0.0259936   0.110111
 -0.0420585  -0.0680521
=#
```
The `deigen` function return an object of type `DEigen{Float64, 2}`. This obeject can be unpacked to into a tuple of two variables `values_list` and `vectors_list`.
The variable `values_list` stores all $K$ derivatives of all $M$ eigenvalues in the way in a $M\times(K+1)$ array
```math
\begin{bmatrix}
\lambda^{(0)}_{1} & \lambda^{(1)}_{1} & \cdots & \lambda^{(K)}_{1}\\
\lambda^{(0)}_{2} & \lambda^{(1)}_{2} & \cdots & \lambda^{(K)}_{2} \\
\vdots            & \vdots            & \ddots & \vdots  \\
\lambda^{(0)}_{M} & \lambda^{(1)}_{M} & \cdots & \lambda^{(K)}_{M}
\end{bmatrix}
```
and the variable `vectors_list` stores all $K$ derivatives of all $M$ eigenvectors in a $M\times M\times (K+1)$ array
```math
\overset{ \text{vectors\_list[:, :, 1]} }{
  \begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(0)}_1   & v^{(0)}_2   & \cdots & v^{(0)}_M   \\
  \vert & \vert &        & \vert
  \end{bmatrix}
},\,
\overset{ \text{vectors\_list[:, :, 2]} }{
\begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(1)}_1   & v^{(1)}_2   & \cdots & v^{(1)}_M   \\
  \vert & \vert &        & \vert
\end{bmatrix}
}
,
\ldots
,\,
\overset{ \text{vectors\_list[:, :, K+1]} }{
  \begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(K)}_1   & v^{(K)}_2   & \cdots & v^{(K)}_M   \\
  \vert & \vert &        & \vert
\end{bmatrix}.
}
```


### Derivatives of selected eigenvalues and eigenvectors
One can also choose to compute the derivatives of only selected eigenvalues and eigenvectors.
```julia
using DEigen
# A0 = A(1),  A1 = A'(1), A2 = A''(1)
A0 = [1.0 1.0;
      1.0 2.0]

A1 = [1.0 0.0;
      0.0 0.0]

A2 = [0.0 0.0;
      0.0 0.0]

values0, vectors0 = eigen(A0)

# the indicies for the selected eigenvalues and eigenvectors
selected_ind = 1

values_list, vectors_list = deign(values0[selected_ind], vectors0[:,selected_ind], A0, A1, A2)

#=
DEigen{Float64, 2}
values_list:
1×3 Matrix{Float64}:
 0.381966  0.723607  -0.178885

vectors_list:
2×1×3 Array{Float64, 3}:
[:, :, 1] =
 -0.8506508083520399
  0.5257311121191335

[:, :, 2] =
 0.10514622242382668
 0.17013016167040795

[:, :, 3] =
 -0.025993575698632487
 -0.042058488969530655
=#
```
In this case, we select $P$ indices $[i_1,i_2,\ldots,i_P] \subseteq [1,2,\ldots,M]$, and compute the derivatives of the eigenvalues and eigenvectors corresponding to these indices.
The variable `values_list` stores all $K$ derivatives of selected $P$ eigenvalues in the way in a $P\times(K+1)$ array
```math
\begin{bmatrix}
\lambda^{(0)}_{i_1} & \lambda^{(1)}_{i_1} & \cdots & \lambda^{(K)}_{i_1}\\
\lambda^{(0)}_{i_2} & \lambda^{(1)}_{i_2} & \cdots & \lambda^{(K)}_{i_2} \\
\vdots            & \vdots            & \ddots & \vdots  \\
\lambda^{(0)}_{i_P} & \lambda^{(1)}_{i_P} & \cdots & \lambda^{(K)}_{i_P}
\end{bmatrix}
```
and the variable `vectors_list` stores all $K$ derivatives of selected $P$ eigenvectors in a $M\times P\times (K+1)$ array
```math
\overset{ \text{vectors\_list[:, :, 1]} }{
  \begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(0)}_{i_1}   & v^{(0)}_{i_2}   & \cdots & v^{(0)}_{i_P}   \\
  \vert & \vert &        & \vert
  \end{bmatrix}
},\,
\overset{ \text{vectors\_list[:, :, 2]} }{
\begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(1)}_{i_1}   & v^{(1)}_{i_2}   & \cdots & v^{(1)}_{i_P}   \\
  \vert & \vert &        & \vert
\end{bmatrix}
}
,
\ldots
,\,
\overset{ \text{vectors\_list[:, :, K+1]} }{
  \begin{bmatrix}
  \vert & \vert &        & \vert \\
  v^{(K)}_{i_1}   & v^{(K)}_{i_2}   & \cdots & v^{(K)}_{i_P}   \\
  \vert & \vert &        & \vert
\end{bmatrix}.
}
```

## Reference
T. Mach and M. A. Freitag, Solving the Parametric Eigenvalue Problem by Taylor Series and Chebyshev Expansion, preprint (2023)
(https://arxiv.org/pdf/2302.03661)
