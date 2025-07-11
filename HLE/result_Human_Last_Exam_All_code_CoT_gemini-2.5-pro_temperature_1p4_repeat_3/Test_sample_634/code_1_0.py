import numpy as np

# The problem is to find the greatest common right divisor (GCRD) of two polynomial matrices.
# Based on the step-by-step derivation using the Euclidean algorithm (via row reduction),
# the GCRD is the 2x2 identity matrix.

# The final equation is GCRD = D, where D is the resulting matrix.
# Here we define and print the matrix D.
GCRD = np.array([
    [1, 0],
    [0, 1]
])

print("The greatest common right divisor matrix is:")
# We output each number in the final resulting matrix
print(f"[{GCRD[0, 0]}, {GCRD[0, 1]}]")
print(f"[{GCRD[1, 0]}, {GCRD[1, 1]}]")