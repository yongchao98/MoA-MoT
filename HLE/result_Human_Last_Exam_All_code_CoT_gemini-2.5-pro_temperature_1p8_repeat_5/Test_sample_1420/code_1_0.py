import sys
from fractions import Fraction

# --- Parameters ---
# The user specified spatial dimension.
d = 3
# The upper critical dimension for the O(n) universality class (e.g., Ising, XY, Heisenberg models).
d_c = 4
# The number of components of the order parameter.
# We assume the Ising model, which is the most common case where n=1.
n = 1

# --- Calculation ---
# 1. Calculate the epsilon parameter.
epsilon = d_c - d

# 2. Calculate the numerator and denominator of the formula for alpha.
# The formula to first order in epsilon is: α = (4 - n) / (2 * (n + 8)) * ϵ
numerator = 4 - n
denominator = 2 * (n + 8)

# 3. Calculate alpha.
alpha = (numerator / denominator) * epsilon

# --- Output ---
print(f"Finding the specific heat scaling exponent α for a system in d={d} dimensions.")
print("This calculation uses the first-order epsilon (ϵ) expansion near the upper critical dimension.")
print(f"Assuming the Ising model universality class (n={n}), the upper critical dimension is d_c={d_c}.")
print("\nFirst, we calculate ϵ:")
print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}\n")

print("The formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")
print("Substituting the values n=1 and ϵ=1:")
# The user requested to see each number in the final equation.
# The sys.stdout.write and flush commands are used for staged printing to be clear.
sys.stdout.write(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}\n")
sys.stdout.flush()
sys.stdout.write(f"α = {numerator} / (2 * {n+8}) * {epsilon}\n")
sys.stdout.flush()
sys.stdout.write(f"α = {numerator} / {denominator} * {epsilon}\n")
sys.stdout.flush()

# Represent as a fraction for precision
alpha_fraction = Fraction(numerator, denominator) * epsilon
print(f"α = {alpha_fraction.numerator}/{alpha_fraction.denominator}\n")

print(f"The final calculated value for α is approximately: {alpha}")

<<<0.16666666666666666>>>