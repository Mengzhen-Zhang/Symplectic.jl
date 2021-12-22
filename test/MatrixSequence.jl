using Symplectic
a1, a2, a3 = [1 0; 1 1], [2 3; 4 5], [7 8;9 10]
b1, b2, b3, b4 = [6 3; 9 4], [1 7; 3 2], [4 3; 5 2], [6 1; 8 3]
sequence = [a1, a2, a3]
interspersion = [b1, b2, b3, b4]
ms = MatrixSequence( [a1, a2, a3] )
ms2 = MatrixSequence(sequence, interspersion)
replaceSequence(ms2, ms, 2)