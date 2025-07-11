import numpy as np

# Step 1: Define the given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Set up the system of linear equations A * pi = b for the steady state.
# We use the first balance equation, the third, the fourth, and the normalization condition.
# The vector of unknowns is pi = [P0, P1, P2, P3]
#
# Equations:
# 1) -lambda_01 * P0 + lambda_10 * P1 = 0
# 2) lambda_12 * P1 - (lambda_21 + lambda_23) * P2 = 0
# 3) lambda_23 * P2 - lambda_31 * P3 = 0
# 4) P0 + P1 + P2 + P3 = 1

A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Step 3: Solve the system for the probabilities pi
try:
    pi = np.linalg.solve(A, b)
    P0_inf = pi[0]
    P1_inf = pi[1]

    # Step 4: Calculate the required sum
    result = P0_inf + P1_inf

    # Step 5: Print the final result, showing each number in the equation
    print(f"The steady-state probabilities are:")
    print(f"P0(inf) = {P0_inf}")
    print(f"P1(inf) = {P1_inf}")
    print(f"P2(inf) = {pi[2]}")
    print(f"P3(inf) = {pi[3]}")
    print("\nThe required sum is P0(inf) + P1(inf):")
    print(f"{P0_inf} + {P1_inf} = {result}")
    
    # Final answer in the required format
    print(f"\n<<<{result}>>>")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
