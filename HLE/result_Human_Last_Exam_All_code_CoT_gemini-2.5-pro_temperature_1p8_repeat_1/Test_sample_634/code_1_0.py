import numpy as np

# The problem is to find the greatest common right divisor (GCRD)
# of two polynomial matrices P1 and P2.
# Through the matrix Euclidean algorithm (performing row operations on the stacked matrix [P1; P2]),
# we found that the GCRD is the 2x2 identity matrix.

# The elements of the GCRD matrix G(s)
g11 = 1
g12 = 0
g21 = 0
g22 = 1

# Create the matrix using a list of lists
G = [[g11, g12],
     [g21, g22]]

# Print the final result
print("The greatest common right divisor is the matrix G(s) where:")
print(f"g11(s) = {G[0][0]}")
print(f"g12(s) = {G[0][1]}")
print(f"g21(s) = {G[1][0]}")
print(f"g22(s) = {G[1][1]}")
print("\nSo the GCRD matrix is:")
# Use numpy for pretty printing the matrix, though it's not strictly necessary.
# For simple cases, manual formatting is also fine.
try:
    print(np.array(G))
except ImportError:
    print(f"[[{g11}, {g12}],")
    print(f" [{g21}, {g22}]]")
