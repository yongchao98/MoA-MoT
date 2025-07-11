import math

# Coefficients derived from the analysis
A = 9.0
B = -9.0
C = 1.0
exponent_val = 1.0/3.0

# Rounding the coefficients to two decimal places for the output expression
A_r = round(A, 2)
B_r = round(B, 2)
C_r = round(C, 2)
exponent_r = round(exponent_val, 2)

# Print the final analytical expression, outputting each number
print("The approximate analytical expression for y(x) in the large x regime is:")
print(f"y(x) = {A_r:.2f} + {B_r:.2f} * ({C_r:.2f} - x)**{exponent_r:.2f}")
