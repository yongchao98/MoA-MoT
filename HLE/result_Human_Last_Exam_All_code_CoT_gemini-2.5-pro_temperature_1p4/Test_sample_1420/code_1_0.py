import sympy

# Step 1: Define the parameters for the problem.
d = 3  # Spatial dimension
d_c = 4  # Upper critical dimension
n = 1  # Number of order parameter components (assuming Ising universality class)

# Step 2: Calculate epsilon (ϵ)
epsilon = d_c - d

# Step 3: Use the first-order epsilon expansion formula for the specific heat exponent α.
# The formula is α = (4 - n) / (2 * (n + 8)) * ϵ
numerator = 4 - n
denominator = 2 * (n + 8)

# The result as a fraction for precision
alpha_fraction = sympy.Rational(numerator, denominator) * epsilon
# The result as a floating-point number
alpha_float = float(alpha_fraction)

# Step 4: Print the explanation, the formula, the calculation, and the result.
print("Calculation of the specific heat scaling exponent α using first-order ϵ-expansion:")
print("-" * 70)
print(f"Spatial dimension, d = {d}")
print(f"Upper critical dimension, d_c = {d_c}")
print(f"Number of order parameter components, n = {n} (Ising universality class)")
print(f"\nFirst, we calculate epsilon: ϵ = d_c - d = {d_c} - {d} = {epsilon}")

print("\nThe first-order formula for α is: α = (4 - n) / (2 * (n + 8)) * ϵ")

print("\nSubstituting the values into the formula:")
# The prompt requires outputting each number in the final equation.
print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
print(f"α = {numerator} / ({2} * {n+8}) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")

print(f"\nThe result is: α = {alpha_fraction}")
print(f"As a decimal, this is approximately: α ≈ {alpha_float:.4f}")

# The final answer is the floating point value.
final_answer = alpha_float