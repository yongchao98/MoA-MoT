import numpy as np

# The problem asks for the eigenvalue of the neutralino mass matrix
# that is not proportional to the adjustable parameters M_1, M_2, or mu.
# As derived in the plan, this scenario leads to a 2x2 sub-matrix for the
# tilde-Z and tilde-H_a states. To find an eigenvalue independent of M_1 and mu,
# we consider the limit where M_1 and mu are negligible compared to M_Z.

# The mass of the Z boson in GeV.
M_Z = 91.1876

# In the limit M_1 -> 0 and mu -> 0, the sub-matrix becomes:
M_sub = np.array([[0, M_Z],
                  [M_Z, 0]])

# The eigenvalues of this matrix are found by solving the characteristic equation:
# det(M_sub - lambda*I) = 0
# (-lambda)(-lambda) - M_Z*M_Z = 0
# lambda^2 - M_Z^2 = 0
print("The characteristic equation for the relevant sub-matrix simplifies to: lambda^2 - M_Z^2 = 0")
print(f"Substituting the value of M_Z: lambda^2 - {M_Z:.4f}^2 = 0")
print(f"This is equivalent to: lambda^2 - {M_Z**2:.4f} = 0")

# We can solve this using numpy's eigenvalue calculator.
eigenvalues = np.linalg.eigvals(M_sub)

print(f"\nThe two eigenvalues that are not proportional to M_1 or mu are {eigenvalues[0]:.4f} and {eigenvalues[1]:.4f}.")

# The question asks for "the" eigenvalue, singular. It is conventional to refer
# to the positive value, which represents the mass scale.
final_eigenvalue = max(eigenvalues)

print(f"\nSelecting the positive value, the computed eigenvalue is: {final_eigenvalue:.4f}")