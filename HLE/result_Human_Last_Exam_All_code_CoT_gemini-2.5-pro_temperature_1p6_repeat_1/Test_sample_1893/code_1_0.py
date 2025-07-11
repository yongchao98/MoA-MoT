import numpy as np

# In the scenario described, the problem of finding the neutralino mass eigenvalues
# simplifies significantly. We focus on the part of the mass matrix that does not
# directly yield eigenvalues proportional to M1 or mu. This leads us to consider
# a 2x2 sub-matrix in the limit where M1 and mu are zero.

# The mass of the Z boson in GeV.
M_Z = 91.1876

# The simplified sub-matrix for the Z-ino and Higgsino states is:
# [[0, M_Z],
#  [M_Z, 0]]
matrix_to_diagonalize = np.array([[0.0, M_Z],
                                  [M_Z, 0.0]])

# We calculate the eigenvalues of this matrix.
eigenvalues = np.linalg.eigvals(matrix_to_diagonalize)

# The characteristic equation for this matrix is det([[0-λ, M_Z], [M_Z, 0-λ]]) = 0,
# which simplifies to λ^2 - M_Z^2 = 0.
print(f"The characteristic equation to solve is: lambda^2 - ({M_Z})^2 = 0.")

# The solutions (eigenvalues) are +M_Z and -M_Z.
print(f"The two eigenvalues are {eigenvalues[0]:.4f} and {eigenvalues[1]:.4f}.")

# The question asks for the eigenvalue not proportional to M1, M2, or mu.
# We choose the positive value, corresponding to the mass.
positive_eigenvalue = max(eigenvalues)

print(f"\nThe required eigenvalue is M_Z, which is approximately {positive_eigenvalue:.4f}.")
