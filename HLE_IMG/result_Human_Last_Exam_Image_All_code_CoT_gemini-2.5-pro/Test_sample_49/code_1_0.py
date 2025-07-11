import math

# This script presents the derived formula for the cutoff frequency at node a0.

# The analysis of the circuit shows that the cutoff frequency f_c follows the formula:
# f_c = 1 / (2 * pi * R_th * C)
# where R_th is the Thevenin equivalent resistance.

# The Thevenin resistance R_th was derived as:
# R_th = r * (1 + sqrt(3))

# Substituting R_th into the frequency formula gives the final expression.
# The numbers present in this final equation are 1, 2, 3, and pi.

# Define the symbolic and numeric components of the equation
num_one = 1
num_two = 2
num_three = 3
pi_symbol = "pi"
r_symbol = "r"
C_symbol = "C"

# Print the final equation in a clear, symbolic format.
print("The cutoff frequency (f_c) at node a0 is given by the equation:")
print(f"f_c = {num_one} / ({num_two} * {pi_symbol} * ({num_one} + sqrt({num_three})) * {r_symbol} * {C_symbol})")

# For practical calculation, we can show the numerical value of the constant coefficient.
resistance_factor = num_one + math.sqrt(num_three)
denominator_coefficient = num_two * math.pi * resistance_factor
print("\nNumerically, the equation is approximately:")
print(f"f_c = {num_one} / ({denominator_coefficient:.4f} * {r_symbol} * {C_symbol})")
