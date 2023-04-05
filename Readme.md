
# Table of Contents

1.  [This document is under construction.](#orgbbb13a2)
2.  [Symplectic.jl](#org2322aa1)
    1.  [Mathematical functions](#org2a26651)
        1.  [Symplectic form](#org0bc84de)
        2.  [Basis of the symplectic vector space](#org7cd6e92)
        3.  [Symplectic cayley transform](#orge7c3d99)
        4.  [Pre-Iwasawa factorization](#org642c433)
        5.  [Utilities](#org97f64b2)
    2.  [Physical functions](#org11a2637)
        1.  [Phase shifting](#orgd37b284)
        2.  [Beam-splitter](#org496897e)
        3.  [Two-mode squeezing](#org401929d)
        4.  [Circulator](#org84d6cb0)
        5.  [Adpative control](#org8fc5bf1)
        6.  [Interferometer](#orge086a2f)
        7.  [Dilation](#orga17a4ae)
        8.  [Gaussian channel](#orgf05dcf2)
        9.  [Squeezed vacuum](#orgaa434bc)
        10. [Coupled-mode system](#org1813e45)
    3.  [Circuit](#orgd9fb8f0)


<a id="orgbbb13a2"></a>

# This document is under construction.


<a id="org2322aa1"></a>

# Symplectic.jl

A Julia package for construction and manipulation of symplectic matrices, Gaussian channels (as in quantum information), and circuits comprised of symplectic matrices and Gaussian channels.


<a id="org2a26651"></a>

## Mathematical functions


<a id="org0bc84de"></a>

### Symplectic form

Throughout the document, the skew-symmetric form defining the symplectic matrix is of the following form
$$\Omega = \oplus_{i=1}^n \begin{pmatrix}
0 & 1 \\
-1 & 0
\end{pmatrix},$$
i.e., a 2n-by-2n block-diagonal matrix with n identical sub-blocks. This matrix is provided as a struct `SymplecticForm` which resembles the `UniformScaling` struct in `LinearAlgebra`. Similar as the symbol `I` for `UniformScaling`, the symbol `$\Omega$` is reserved for `SymplecticForm`. One can quickly create the skew-symmetric matrix as shown above by calling `$\Omega$(2n)` where `2n` represents the size of the matrix. Just as `I`, `$\Omega$` can be used without an argument. That is to say, basis matrix arithmetics for `$\Omega$` are supported so long as the size of the matrix can be inferred from the context.


<a id="org7cd6e92"></a>

### Basis of the symplectic vector space

With the explicit form of `$\Omega$` as is specified above, the basis of the underlying symplectic vector space is accordingly fixed by default. We refer to this default choice of symplectic basis the **QPQP** basis. However, for certain use cases, it is more convenient to work with an alternative basis, referred to as the **QQPP** basis, in which the symplectic form is of the following form
$$
\Omega = \begin{pmatrix}
 0_n & I_n \\
 -I_n & 0_n
\end{pmatrix}
$$
with $0_n$ represents an n-by-n zero matrix and $I_n$ an n-by-n identity matrix. The function `toQQPPBasis` (respectively, `toQPQPBasis`) is provided to convert a matrix represented in the QPQP (respectively, QQPP) basis to the QQPP (respectively, QPQP) basis.


<a id="orge7c3d99"></a>

### Symplectic cayley transform

Cayley transform is an important technique to parametrize orthogonal or unitary matrices. With appropriate modification, it can also be used to parametrize other groups including the group of symplectic matrices.

The symplectic cayley transform is defined as follows
$$S = (\Omega M - I/2)(\Omega M + I/2)^{-1}$$
where `M` is a 2n-by-2n **real** matrix and `S` is a 2n-by-2n symplectic matrix; and its inverse is
$$M = \Omega (S - I)^{-1}(S + I)/2$$
i.e., the `M` thus obtained should be the same as that in the equation above.

The symplectic cayley transform is implemented as a single-argument funciton `symplecticCayleyTransform` (with alias `cayley`), and its inverse as `inverseSymplecticCayleyTransform` (with alias `invcaylay`).


<a id="org642c433"></a>

### Pre-Iwasawa factorization

A symplectic matrix $S$ can always be written as the multiplication of three symplectic matrices of specific structures, which is known as the pre-Iwasawa factorization. The factorization can be specificied by the following equation (note that the matrices are represented in the **QQPP** basis) 
\\[
 S = \begin{pmatrix}
    I<sub>n</sub> & 0<sub>n</sub>   
    P & I<sub>n</sub>
 \end{pmatrix}

\begin{pmatrix}
   L & 0_n \\
   0_n & (L^T)^{-1}
\end{pmatrix}

 \begin{pmatrix}
    X & Y   
   -Y & X
 \end{pmatrix}.
\\]
Here, $P$ is an n-by-n symmetric matrix; $L$ is an n-by-n non-singular matrix; $X$ and $Y$ are n-by-n matrices satisfying that $X+iY$ is a unitary matrix, i.e., the right-most of the three matrices on the right-hand-side of the above equation is an orthogonal matrixb. The pre-Iwasawa factorization is unique and well-defined for all symplectic matrices.

The pre-Iwasawa factorization is implemented as the function `preIwasawaFactorization` (or `preiwa`), which yields a tuple of three symplectic matrices as are on the right-hand-side of the above equation.


<a id="org97f64b2"></a>

### Utilities

We provide `hs_norm` to facilliate the calculation of the Hibert-Schmidt norm of a matrix. That is, `hs_norm(M)` yields $\mbox{tr}(M^T*M)$.

To check whether a matrix is symplectic, the user can use the function `nonSymplecticity`. It takes a matrix $S$ as the input and yields the Hilber-Schmidt norm of $S^T \Omega S - \Omega$.

The function `dsum` (or the binary operator `$\oplus$`) is used to produce direct sum of two matrices. It can also take more than two matrices as input: `dsum(A1, A2, ...)`

The binary operator `\otimes` is assigned to `Base.kron` to allow tensor product.


<a id="org11a2637"></a>

## Physical functions

In physics, symplectic matrices significantly simplifies the analysis of unitary Gaussian processes transforming quadratures operators into their linear combinations. Representations of some commonly used physical componets are provided.


<a id="orgd37b284"></a>

### Phase shifting

The 2-by-2 symplectic matrix corresponding to a single-mode phase-shifting operation can be created by the `phaseShifting` with a single real argument. When multiple arguments are given, e.g. ~phaseShifting(x1, x2, x3), it yields the symplectic matrix correpsonding to the simultaneous application of multiple single-mode operations with angles specificed correspondingly by the 
arguments.


<a id="org496897e"></a>

### Beam-splitter

The function `beamSplitter` yields a 4-by-4 symplectic matrix corresponding to a two-mode beam-splitter, with the angle specificed by the input argument.

The alternative method `beamSplitter(angle, m1, m2, n)` yields a beam-splitter symplectic matrix between mode m1 and m2 in an n-mode system. When n is omitted, the `max(m1, m2)` will be used as the total number of modes in the system.


<a id="org401929d"></a>

### Two-mode squeezing

The symplectic matrix representation of a two-mode squeezing operation can be created by the function `amplifier` with its argument specifing the gain coefficient.

When called with `amplifier(G, m1, m2, n)`, it yields the two-mode squeezing operation between the mode m1 and the mode m2 in a n-mode system. When `n` is omitted, `max(m1, m2)` will be used as the number of modes of the whole system.


<a id="org84d6cb0"></a>

### Circulator

The function `cirulator(perm::Vector)` yields a symplectic matrix of a circulator, i.e., a permutation of modes specified by the vector `perm`. Alternatively, it supports the method `cirulator(perm...)`.


<a id="org8fc5bf1"></a>

### Adpative control

The function `teleport(S::AbstractMatrix, inModes::Vector, outModes::Vector)` implements of the main result in Phys. Rev. Lett. 120, 020502<sup><a id="fnr.1" class="footref" href="#fn.1" role="doc-backlink">1</a></sup>. Here, `S` represents a 2n-by-2n symplectic matrix representing a unitary Gaussian operatio on an n-mode system. `inModes` and `outModes` are vectors of equal sizes representing the input and output modes. Let `ancModes` denote the vector of the modes that are not in `inModes`, and `idlModes` the vector of the modes that are not in `outModes`. Let `Usq` denote `2*ancModes.-1`, and `Hm` denote `2*idlModes.-1`, `In` denote the set of indices either in `2*inModes.-1` or `2*inModes`, and `Out` denote the set of indices either in `2*outModes.-1` or `2*outModes`. This function outputs `S[Out,In]-S[Out,Usq]*(S[Hm,Usq])^(-1)*S[Hm,In]` which is guaranteed to be a symplectic matrix so long as `S[Hm,Usq]` is non-singular.

The function `feedforward(S::AbstractMatrix, inModes::Vector, outModes::Vector)` yields the product `-S[Out,USq]*(S[Hm,Usq])^(-1)` directly.

The function `adaptiveMeasurement(F::AbstractMatrix, outModes::Vector, n::Integer)` yields the matrix defined in Eq.(4.3.11) of <sup><a id="fnr.2" class="footref" href="#fn.2" role="doc-backlink">2</a></sup> in the default QPQP basis. Here `F` is a linear map from `Hm` to `outModes`.


<a id="orge086a2f"></a>

### Interferometer

The function `interferenceBasedSequence(S::AbstractMatrix; T=I(4))` yields an array of symplectic matrices consisting of multiple copies of `S` interspersed with symplectic matrices that are direct sum of single-mode diagonal blocks. The product of the array is equal to `T`.

The method `interferenceBasedSequence(Ss; T=I(4))` allows replacing single `S` with a sequnce of symplectic matrices.

This function has an alias `infseq`. It implements the main reult in <sup><a id="fnr.3" class="footref" href="#fn.3" role="doc-backlink">3</a></sup>.


<a id="orga17a4ae"></a>

### Dilation

The function `dilate(S::AbstractMatrix)` implements Theorem 3.15 in <sup><a id="fnr.2.100" class="footref" href="#fn.2" role="doc-backlink">2</a></sup>. The input can be an arbitrary real square matrix; and the output is a symplectic matrix.


<a id="orgf05dcf2"></a>

### Gaussian channel

The function `channel(S::AbstractMatrix, Venv::AbstractMatrix, inModes::Vector, outModes::Vector)` yields a pair of matrices representing a Gaussian channel transforming quantum states of `inModes` to states of `outModes`.

The symplectic matrix `S` represents the unitary Gaussian operation on both the system modes and the environment modes. `Venv` represents the covariance matrix of the environment modes.


<a id="orgaa434bc"></a>

### Squeezed vacuum

The function `squeezedVacuum(x::Vector)` generates the covariance matrix `length(x)` modes with their Q-quadrature squeezed. For the ith mode, the degree of squeezing is `x[i]` dB.


<a id="org1813e45"></a>

### Coupled-mode system


<a id="orgd9fb8f0"></a>

## Circuit


# Footnotes

<sup><a id="fn.1" href="#fnr.1">1</a></sup> <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.020502>

<sup><a id="fn.2" href="#fnr.2">2</a></sup> <https://arxiv.org/pdf/2107.01474v1.pdf>

<sup><a id="fn.3" href="#fnr.3">3</a></sup> <https://www.nature.com/articles/s41534-022-00581-9>
