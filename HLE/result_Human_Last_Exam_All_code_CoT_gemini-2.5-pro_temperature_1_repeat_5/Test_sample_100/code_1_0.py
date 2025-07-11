import numpy as np

# Based on the mathematical derivation, the integral evaluates to the expression:
# (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1

# Define the coefficients and constants from the final expression
c1_num = 8
c1_den = 15
c2_num = 1
c2_den = 3
c3_num = -1
c3_den = 2
c4 = 1

# Calculate the numerical value of the expression
pi = np.pi
result_value = (c1_num / c1_den) * pi**8 + \
               (c2_num / c2_den) * pi**2 + \
               (c3_num / c3_den) * pi + \
               c4

# Output the final equation with each number, as requested.
print("The final equation for the integral is:")
print(f"({c1_num}/{c1_den})*pi^8 + ({c2_num}/{c2_den})*pi^2 + ({c3_num}/{c3_den})*pi + {c4}")
print(f"\nNumerical value: {result_value}")
