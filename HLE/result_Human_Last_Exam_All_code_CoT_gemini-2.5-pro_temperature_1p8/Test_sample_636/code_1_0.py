import numpy as np

# Step 1: Define the system matrices A and B
A = np.array([[-1, 1], [1, 0]])
B = np.array([[1, 2], [1, 0]])

# Step 2: Define the desired eigenvalues and find the desired characteristic polynomial
# Desired eigenvalues are p = -1 +/- j
# Desired characteristic polynomial is (s - (-1+j))(s - (-1-j)) = s^2 + 2s + 2
# From the polynomial, the desired trace is -2 and the desired determinant is 2.
desired_trace = -2
desired_det = 2
print("Pole Placement Problem\n")
print(f"A = \n{A}\n")
print(f"B = \n{B}\n")
print(f"Desired eigenvalues: -1 + 1j, -1 - 1j")
print(f"Desired characteristic polynomial: s^2 + 2s + 2")
print("-" * 30)

# Step 3 & 4: Set up the equations for the elements of F = [[f1, f2], [f3, f4]]
# The characteristic polynomial of the closed-loop system A_cl = A + BF is:
# s^2 - trace(A+BF)s + det(A+BF)
# Equating coefficients with the desired polynomial gives:
# 1) trace(A+BF) = -2
# 2) det(A+BF) = 2
#
# trace(A+BF) = trace(A) + trace(BF) = (-1+0) + trace([[f1+2*f3, f2+2*f4], [f1, f2]])
# -1 + f1 + 2*f3 + f2 = -2  =>  f1 + f2 + 2*f3 = -1  (Equation 1)
#
# This gives a system of 2 equations and 4 unknowns. We can simplify by choosing
# two variables. Let's find a simple solution by setting f3=0 and f4=0.

# Step 5: Solve the simplified system of equations
# With f3=0 and f4=0, F = [[f1, f2], [0, 0]]
# The closed-loop matrix A_cl = A + B @ [[f1, f2], [0, 0]] becomes:
# A_cl = [[-1, 1], [1, 0]] + [[f1, f2], [f1, f2]] = [[-1+f1, 1+f2], [1+f1, f2]]
#
# Now the equations for trace and determinant become:
# trace(A_cl) = (-1+f1) + f2 = f1 + f2 - 1
# det(A_cl) = f2(-1+f1) - (1+f2)(1+f1) = -f2 + f1*f2 - (1 + f1 + f2 + f1*f2) = -f1 - 2*f2 - 1
#
# Set these to the desired values:
# 1) f1 + f2 - 1 = -2  =>  f1 + f2 = -1
# 2) -f1 - 2*f2 - 1 = 2  => -f1 - 2*f2 = 3
#
# We solve this 2x2 system for f1 and f2.
# C * f_vec = d
C = np.array([[1, 1], [-1, -2]])
d = np.array([-1, 3])
f1_f2 = np.linalg.solve(C, d)

f1 = f1_f2[0]
f2 = f1_f2[1]
f3 = 0.0
f4 = 0.0

# The calculated feedback gain matrix F
F = np.array([[f1, f2], [f3, f4]])

print("By setting f3=0 and f4=0, we solved for f1 and f2.")
print(f"The calculated state feedback gain matrix F is:")
# The instruction was "output each number in the final equation"
# We print the matrix F that completes the equation A_cl = A + B@F
print(F)
print("-" * 30)


# Step 6: Verification
print("Verification:")
# Calculate the closed-loop system matrix A_cl
A_cl = A + B @ F

print("The closed loop equation is A_cl = A + B @ F")
print("A_cl = \n", A, "\n+ \n", B, "\n@\n", F)
print("\nA_cl = \n", A_cl)


# Calculate eigenvalues of the closed-loop system
eigenvalues = np.linalg.eigvals(A_cl)
print(f"\nThe eigenvalues of the closed-loop system A+BF are: {np.round(eigenvalues, 5)}")
print("This matches the desired eigenvalues.")
