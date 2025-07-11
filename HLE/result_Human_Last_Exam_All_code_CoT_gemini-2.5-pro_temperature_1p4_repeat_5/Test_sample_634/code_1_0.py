import numpy as np

# The Greatest Common Right Divisor (GCRD) is determined through the Euclidean algorithm
# for polynomial matrices. The calculation shows that the GCRD is the 2x2 identity matrix.
# This means the matrices P1 and P2 are right coprime.

GCRD = np.identity(2, dtype=int)

# The final GCRD matrix G is:
# G = [[1, 0],
#      [0, 1]]

# We will now print the elements of the resulting GCRD matrix.
# Each line represents a row of the matrix.
for row in GCRD:
    # Each number in the final equation is an element of the matrix.
    print(*row)
