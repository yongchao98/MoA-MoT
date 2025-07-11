import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# To find the stationary probabilities P_inf = [pi_0, pi_1, pi_2, pi_3],
# we set the derivatives of the system to zero, which gives a system of linear equations.
# We select three independent equations from the system and add the normalization condition.
#
# 1) -lambda_01*pi_0 + lambda_10*pi_1 = 0
# 2) lambda_12*pi_1 - (lambda_21 + lambda_23)*pi_2 = 0
# 3) lambda_23*pi_2 - lambda_31*pi_3 = 0
# 4) pi_0 + pi_1 + pi_2 + pi_3 = 1
#
# This can be written in matrix form A * P_inf = b

A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations for the stationary probabilities
P_inf = np.linalg.solve(A, b)

# Extract the probabilities for state 0 and state 1
P0_inf = P_inf[0]
P1_inf = P_inf[1]

# Calculate the required sum
result = P0_inf + P1_inf

# As requested, output each number in the final equation.
# The final equation is P0(+infinity) + P1(+infinity) = result
print(f"P0(+infinity) = {P0_inf}")
print(f"P1(+infinity) = {P1_inf}")
print(f"The equation is: {P0_inf} + {P1_inf} = {result}")
