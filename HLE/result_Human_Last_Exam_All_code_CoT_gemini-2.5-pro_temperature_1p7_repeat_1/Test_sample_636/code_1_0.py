import numpy as np

# This script calculates a state feedback gain F to place the poles of A + B*F.
#
# The characteristic equation for the simplified closed-loop system A_cl = A + B*F,
# with F = [[f1, f2], [0, 0]], is derived as follows:
#
# B*F = [[1, 2],   *  [[f1, f2],  =  [[f1, f2],
#        [1, 0]]      [0, 0]]       [f1, f2]]
#
# A_cl = A + B*F = [[-1, 1], + [[f1, f2], = [[f1 - 1, f2 + 1],
#                    [ 1, 0]]   [f1, f2]]   [f1 + 1,      f2]]
#
# The characteristic polynomial is det(λI - A_cl):
# |λ - (f1-1)   -(f2+1) |
# |-(f1+1)     λ - f2   |  = (λ - f1 + 1)(λ - f2) - (f1 + 1)(f2 + 1)
#                          = λ^2 + λ*(1 - f1 - f2) - (1 + f1 + 2*f2)
#
# The desired polynomial from eigenvalues -1 ± j is λ^2 + 2λ + 2.
#
# Equating coefficients:
# λ^1: 1 - f1 - f2 = 2        => f1 + f2 = -1
# λ^0: -(1 + f1 + 2*f2) = 2    => f1 + 2*f2 = -3
#
# We solve this 2x2 system for f1 and f2.
# Let C be the coefficient matrix and d be the constant vector.
# [[1, 1], [1, 2]] * [[f1], [f2]] = [[-1], [-3]]
C = np.array([[1., 1.], [1., 2.]])
d = np.array([-1., -3.])

# Solve for f = [f1, f2]
f_params = np.linalg.solve(C, d)
f1 = f_params[0]
f2 = f_params[1]

# Construct the full gain matrix F
F = np.array([[f1, f2],
              [0., 0.]])

# Define system matrices for verification
A = np.array([[-1., 1.],
              [1., 0.]])
B = np.array([[1., 2.],
              [1., 0.]])

# Calculate the closed-loop system matrix A_cl = A + B*F
A_cl = A + B @ F

# Verify the eigenvalues of the closed-loop system
eigenvalues = np.linalg.eigvals(A_cl)

print("This script finds one possible gain matrix F by assuming the second control input is not used (f3=f4=0).\n")

print("The state matrix A is:")
print(A)
print("\nThe input matrix B is:")
print(B)
print("\nThe desired eigenvalues are: -1 + 1j and -1 - 1j\n")

print("The calculated state feedback gain F is:")
print(F)

print("\nThe final equation for the closed-loop matrix is A_cl = A + B*F.")
print(f"\n{A}\n\n+\n\n{B}\n\n*\n\n{F}\n\n=\n\n{A_cl}\n")


print("\nTo verify, the eigenvalues of the resulting matrix A_cl are:")
# Using np.round to clean up display of floating point numbers
print(np.round(eigenvalues, decimals=5))