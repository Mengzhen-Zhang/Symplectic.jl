
# Table of Contents

1.  [This document is under construction.](#orgbe30199)
2.  [Symplectic.jl](#org1ef4842)
    1.  [Mathematical functions](#orgc2c34e4)
        1.  [Symplectic form](#orga837937)
        2.  [Symplectic cayley transform](#org52a247d)


<a id="orgbe30199"></a>

# This document is under construction.


<a id="org1ef4842"></a>

# Symplectic.jl

A Julia package for construction and manipulation of symplectic matrices, Gaussian channels (as in quantum information), and circuits comprised of symplectic matrices and Gaussian channels.


<a id="orgc2c34e4"></a>

## Mathematical functions


<a id="orga837937"></a>

### Symplectic form

Throughout the document, the skew-symmetric form defining the symplectic matrix is of the following form
$$\Omega = \oplus_{i=1}^n \begin{pmatrix}
0 & 1 \\
-1 & 0
\end{pmatrix},$$
i.e., a 2n-by-2n block-diagonal matrix with n identical sub-blocks. This matrix is provided as a struct `SymplecticForm` which resembles the `UniformScaling` struct in `LinearAlgebra`. Similar as the symbol `I` for `UniformScaling`, the symbol `$\Omega$` is reserved for `SymplecticForm`. One can quickly create the skew-symmetric matrix as shown above by calling `$\Omega$(2n)` where `2n` represents the size of the matrix. Just as `I`, `$\Omega$` can be used without an argument. That is to say, basis matrix arithmetics for `$\Omega$` are supported so long as the size of the matrix can be inferred from the context.


<a id="org52a247d"></a>

### Symplectic cayley transform

Cayley transform is an important technique to parametrize orthogonal or unitary matrices. With appropriate modification, it can also be used to parametrize other groups including the group of symplectic matrices.

The symplectic cayley transform is defined as follows
$$S = (\Omega M - I/2)(\Omega M + I/2)^{-1}$$
where `M` is a 2n-by-2n **real** matrix and `S` is a 2n-by-2n symplectic matrix; and its inverse is
$$M = \Omega (S - I)^{-1}(S + I)/2$$
i.e., the `M` thus obtained should be the same as that in the equation above.

The symplectic cayley transform is implemented as a single-argument funciton `symplecticCayleyTransform` (with alias `cayley`), and its inverse as `inverseSymplecticCayleyTransform` (with alias `invcaylay`).

