
# Table of Contents

1.  [This document is under construction.](#org71300e8)
2.  [Symplectic.jl](#orge929fe6)
    1.  [Mathematical functions](#org1e21683)
        1.  [Symplectic form](#org30c4bab)
        2.  [Symplectic cayley transform](#org407455b)


<a id="org71300e8"></a>

# This document is under construction.


<a id="orge929fe6"></a>

# Symplectic.jl

A Julia package for construction and manipulation of symplectic matrices, Gaussian channels (as in quantum information), and circuits comprised of symplectic matrices and Gaussian channels.


<a id="org1e21683"></a>

## Mathematical functions


<a id="org30c4bab"></a>

### Symplectic form

Throughout the document, the skew-symmetric form defining the symplectic matrix is of the following form
$$\Omega = \oplus_{i=1}^n \begin{pmatrix}
0 & 1 \\
-1 & 0
\end{pmatrix},$$
i.e., a 2n-by-2n block-diagonal matrix with n identical sub-blocks. This matrix is provided as a struct `SymplecticForm` which resembles the `UniformScaling` struct in `LinearAlgebra`. Similar as the symbol `I` for `UniformScaling`, the symbol `\Omega` is reserved for `SymplecticForm`. One can quickly create the skew-symmetric matrix as shown above by calling `\Omega(2n)` where `2n` represents the size of the matrix. Just as `I`, `\Omega` can be used without an argument. That is to say, basis matrix arithmetics for `\Omega` are supported so long as the size of the matrix can be inferred from the context.


<a id="org407455b"></a>

### Symplectic cayley transform

Cayley transform is an important technique to parametrize orthogonal or unitary matrices. With appropriate modification, it can also be used to parametrize other groups including the group of symplectic matrices.

The symplectic cayley transform is defined as follows
$$S = (\Omega M - I/2)(\Omega M + I/2)^{-1}$$
where `M` is a 2n-by-2n **real** matrix and `S` is a 2n-by-2n symplectic matrix; and its inverse is
$$M = \Omega (S - I)^{-1}(S + I)/2$$
i.e., the `M` thus obtained should be the same as that in the equation above.

The symplectic cayley transform is implemented as a single-argument funciton `symplecticCayleyTransform` (with alias `cayley`), and its inverse as `inverseSymplecticCayleyTransform` (with alias `invcaylay`).

