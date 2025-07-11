import numpy as np

# Based on the problem description, the condition for dynamic enhancement
# requires M1 = M2 and cos(2*beta) = 0.
# We are asked for an eigenvalue not proportional to M1, M2, or mu.
# This points to a special case where we set the adjustable mass parameters to zero:
# M1 = M2 = 0 and mu = 0.
# Under these conditions, the original 4x4 neutralino mass matrix simplifies to:
#
# M_N = [[0, 0,   0,    0],
#        [0, 0,   M_Z,  0],
#        [0, M_Z, 0,    0],
#        [0, 0,   0,    0]]
#
# Two eigenvalues are clearly 0. The other two come from the 2x2 submatrix
# mixing the tilde(Z) and tilde(H_a) states.

# Mass of the Z boson in GeV.
mz_val = 91.1876

# The relevant 2x2 submatrix.
submatrix = np.array([[0, mz_val],
                      [mz_val, 0]])

# The characteristic equation for this submatrix is det(submatrix - lambda*I) = 0.
# This gives: (0 - lambda) * (0 - lambda) - (M_Z) * (M_Z) = 0
# Which simplifies to: lambda^2 - M_Z^2 = 0.

print("The characteristic equation for the relevant eigenvalues (lambda) is:")
# We print each number required for the final equation.
# Equation form: lambda^2 - C = 0
C = mz_val**2
print(f"(lambda)^2 - ({mz_val})^2 = 0")
print(f"Which evaluates to: (lambda)^2 - {C} = 0")

# The solutions are lambda = +M_Z and lambda = -M_Z.
# The question asks to compute the eigenvalue. We provide the positive value.
eigenvalues = np.linalg.eigvals(submatrix)
positive_eigenvalue = max(eigenvalues)

print("\nThe positive eigenvalue that satisfies the condition is:")
print(positive_eigenvalue)