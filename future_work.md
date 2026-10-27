## Some ideas on the new symplectic library on Julia

We seek to combine four modules into the new `Symplectic.jl` library: the original
`Symplectic.jl` library, the `SymplecticDecompositions.jl` module, the
`torontonian_julia` module and the `trace_distance_lanczos` module. Here is a
brief description of the contents of each module.

*  `Symplectic.jl`. Currently has some matrix decompositions, like the pre-Iwasawa
decomposition, and many quantum optics definitions of phase space related objects,
such as the symplectic representation quantum states and channels. There are very
useful matrix functions such as the Cayley transform. There is a lot of effort in
defining the symplectic form in a similar way as the identity matrix is defined in
Julia.
* `SymplecticDecompositions.jl`. It has the decompositions in [Martin’s](https://cdnsciencepub.com/doi/10.1139/cjp-2024-0070) paper.
Documentation is needed.
* `trace_distance_lanczos`. Contains some useful matrix decompositions and
the definition of Bargmann invariants, the Lanczos algorithm implementation of the
computation of the trace distance, and some analytical bounds on the trace distance
between Gaussian states.
* `torontonian_julia`. Contains a `DoubleFloats` implementation of the
torontonian function using both the recursive algorithm and the direct definition. Also
has the function for the computation of threshold detection probabilities of Gaussian states.

Taking into account the previous content, we would like to propose the following structure for
the new `Symplectic.jl` library:

* `SymplecticDecompositions` module. This is an extension of
`SymplecticDecompositions.jl` that includes everything new included in the
Lanczos algorithm project, the Cayley transform of `Symplectic.jl` and the
recursive matrix decomposition used in the computation of the torontonian.
Everything can be coded so that they allow for arbitrary floating point precision (using
`BigFloats`).
* `SpecialMatrixFunctions` module. This module includes code for computing the
(loop) torontonian and the (loop) hafnian (perhaps also the permanent, although the
`BosonSampling.jl` library seems to already have that implementation) with
arbitrary floating point precision. Here we could include techniques for the computation of the hfanian using the finite-difference sieve.
* `GaussianStates` module. This module contains all the phase space formalism that
was in `Symplectic.jl` (states, unitary transformations, channels, etc)
complemented with Bargmann invariants. Extend to linear superpositions of
Gaussian states where possible. Computation of Fock basis elements or
photon-number probabilities may be included here, since it might be useful for people
looking for code that supports arbitrary precision.
* `TraceDistanceGaussian` module. This module will include everything in the
Lanczos algorithm project. Some functions may be combined in order to remove
redundancy.
* `TraceDistanceGaussianBounds` module. Bounds on the trace distance between
Gaussian states. We can investigate if some of these bounds can be extended to
linear superpositions of Gaussian states.