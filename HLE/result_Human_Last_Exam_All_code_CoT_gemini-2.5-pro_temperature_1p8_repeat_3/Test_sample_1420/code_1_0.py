import sys
from fractions import Fraction

# This script calculates the specific heat exponent alpha using the first-order epsilon expansion.

# --- Step 1: Define the parameters ---
# Spatial dimension of the system
d = 3
# Upper critical dimension for the n-vector model
d_c = 4
# Number of components of the order parameter. We assume the Ising model (n=1).
n = 1

# --- Step 2: Calculate epsilon ---
epsilon = d_c - d

# --- Step 3: Use the one-loop formula for alpha ---
# The formula is α = (4 - n) / (2 * (n + 8)) * ε
numerator = 4 - n
denominator = 2 * (n + 8)

alpha_value = (numerator / denominator) * epsilon

# --- Step 4: Print the calculation and result ---
print(f"Calculation of the specific heat exponent α for d={d} using the ε-expansion.")
print("The first-order formula for the n-vector model is: α = (4 - n) / (2 * (n + 8)) * ε")
print(f"Assuming the Ising model (n={n}) and with ε = d_c - d = {d_c} - {d} = {epsilon}, we substitute the values.")
print("\nFinal Equation:")
# The prompt requires showing the equation with each number.
print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
print(f"α = {numerator} / {denominator}")

# Representing the result as a fraction and a float
final_fraction = Fraction(numerator, denominator)
print(f"α = {final_fraction.numerator}/{final_fraction.denominator}")
print(f"α ≈ {alpha_value:.4f}")

# The final answer in the required format for the system.
# The system will extract this line.
sys.stdout.write(f"\n<<<{alpha_value:.4f}>>>")