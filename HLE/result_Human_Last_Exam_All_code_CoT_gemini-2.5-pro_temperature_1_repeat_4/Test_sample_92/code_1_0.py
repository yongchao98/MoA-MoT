import numpy as np

# The problem of finding the escape probability p = p_0 reduces to solving
# a system of two linear equations with two variables: p and p_2023.
# The equations are derived from the master equation for the random walk
# at states n=2022 and n=2023.
#
# Let p be the escape probability for n <= 2022.
# Let p_2023 be the escape probability for n = 2023.
#
# The system of equations is:
# 1) 82 * p - 52 * p_2023 = 9
# 2) 26 * p - 56 * p_2023 = -9

# We represent this system in matrix form Ax = b
A = np.array([
    [82, -52],
    [26, -56]
])

b = np.array([9, -9])

# Print the equations
print("Solving the system of linear equations Ax = b, where x = [p, p_2023]:")
print(f"Equation 1: {A[0,0]} * p + ({A[0,1]}) * p_2023 = {b[0]}")
print(f"Equation 2: {A[1,0]} * p + ({A[1,1]}) * p_2023 = {b[1]}")
print("-" * 20)

# Solve for x = [p, p_2023]
try:
    solution = np.linalg.solve(A, b)
    p_escape = solution[0]
    p_2023 = solution[1]

    # The probability we are looking for is p = p_0
    print(f"The calculated value for p is: {p_escape}")
    print(f"The calculated value for p_2023 is: {p_2023}")
    print("-" * 20)
    print("The probability that the marble escapes is the value of p.")
    print(f"Final Answer: {p_escape}")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")
