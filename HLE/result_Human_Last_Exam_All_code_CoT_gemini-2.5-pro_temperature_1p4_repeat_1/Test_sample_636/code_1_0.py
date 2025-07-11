import numpy as np

# Set print options for better readability
np.set_printoptions(precision=4, suppress=True)

# 1. Define system matrices
A = np.array([[-1, 1],
              [1, 0]])
B = np.array([[1, 2],
              [1, 0]])

# 2. Define desired eigenvalues and characteristic polynomial coefficients
# Desired polynomial is s^2 + 2s + 2
# From s^2 - trace(A_cl)s + det(A_cl), we need:
# trace(A_cl) = -2
# det(A_cl) = 2

# 3. Set up and solve the system of equations for F
# The equations for the elements of F = [[f1, f2], [f3, f4]] are:
# trace: -1 + f1 + f2 + 2*f3 = -2  =>  f1 + f2 + 2*f3 = -1
# det: -1 - f1 - 2*f2 + 2*f2*f3 - 2*f4 - 2*f1*f4 = 2

# To find a simple solution, we assume f3=0 and f4=0.
# The system simplifies to:
# 1) f1 + f2 = -1
# 2) -1 - f1 - 2*f2 = 2  => f1 + 2*f2 = -3

# We can solve this system [[1, 1], [1, 2]] * [f1, f2]^T = [-1, -3]^T
coeffs = np.array([[1, 1], [1, 2]])
vals = np.array([-1, -3])
solution = np.linalg.solve(coeffs, vals)
f1, f2 = solution[0], solution[1]
f3, f4 = 0.0, 0.0 # Our assumption

# 4. Construct the state feedback gain matrix F
F = np.array([[f1, f2],
              [f3, f4]])

print("The calculated state feedback gain matrix F is:")
print(F)

# 5. Verify the solution by calculating the closed-loop system and its eigenvalues
A_cl = A + B @ F
eigenvalues, _ = np.linalg.eig(A_cl)

print("\n--- Verification ---")
print("The final equation A + B * F = A_cl is:\n")
print(f"{A}\n\n+\n\n{B}\n\n*\n\n{F}\n\n=\n\n{A_cl}")
print("\nDesired eigenvalues: [-1.+1.j -1.-1.j]")
print("Resulting eigenvalues of A + BF:", eigenvalues)
