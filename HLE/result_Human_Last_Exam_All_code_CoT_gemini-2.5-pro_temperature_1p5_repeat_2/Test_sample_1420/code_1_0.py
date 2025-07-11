import sys

# Step 1: Define the parameters for the epsilon expansion.
# The upper critical dimension (d_c) for the n-vector model is 4.
d_c = 4
# The spatial dimension (d) of the system is 3.
d = 3
# We will use the Ising model as the canonical example, where the order
# parameter has n=1 component (spin up or spin down).
n = 1

print(f"Calculating the scaling exponent α for the specific heat in d={d} dimensions.")
print("This will be done using the epsilon expansion to first order.")
print("-" * 60)

# Step 2: Calculate the expansion parameter, epsilon.
print(f"The expansion parameter is ε = d_c - d, where d_c = {d_c} is the upper critical dimension.")
epsilon = d_c - d
print(f"ε = {d_c} - {d} = {epsilon}")
print("-" * 60)

# Step 3: Use the first-order epsilon expansion formula for alpha.
print("The first-order formula for α for the n-vector model is:")
print("α = (4 - n) / (2 * (n + 8)) * ε")
print(f"\nFor the Ising model (n={n}), we substitute the values into the formula:")

# Show the substitution of each number into the equation
print(f"α ≈ ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")

# Step 4: Calculate the intermediate and final values.
numerator = 4 - n
denominator_intermediate = n + 8
denominator_final = 2 * denominator_intermediate

print(f"α ≈ {numerator} / (2 * {denominator_intermediate}) * {epsilon}")
print(f"α ≈ {numerator} / {denominator_final} * {epsilon}")

# Calculate the final result for alpha
alpha = numerator / denominator_final * epsilon

# Use fractions module for a more precise representation of 1/6 if needed
from fractions import Fraction
alpha_fraction = Fraction(numerator, denominator_final) * epsilon

print(f"α ≈ {alpha_fraction.numerator}/{alpha_fraction.denominator}")
print(f"As a decimal, the scaling exponent is α ≈ {alpha:.4f}")

# The final answer in the required format will be based on this calculation
# For this script's purpose, we'll store it to be presented at the end.
final_answer = alpha

# This part is just to show how to format the final answer for the platform
# In a real scenario, you wouldn't print the '<<<' characters like this.
# print(f"\nFinal answer in required format: <<<{final_answer}>>>", file=sys.stderr)
