import sympy

# Define the variables
n = sympy.Integer(1)  # Number of components of the order parameter (n=1 for Ising model)
d = sympy.Integer(3)  # Spatial dimension
d_c = sympy.Integer(4)  # Upper critical dimension

# --- Step 1: Calculate epsilon ---
print(f"The epsilon parameter (ϵ) is defined as d_c - d.")
epsilon = d_c - d
print(f"ϵ = {d_c} - {d} = {epsilon}")
print("-" * 20)

# --- Step 2: Calculate the scaling exponent alpha (α) ---
print(f"The first-order epsilon expansion for the specific heat exponent α is:")
print(f"α = (4 - n) / (2 * (n + 8)) * ϵ")
print("-" * 20)

print(f"Substituting the values n = {n}, and ϵ = {epsilon}:")

# Perform the calculation using symbolic math for a precise fraction
numerator = 4 - n
denominator = 2 * (n + 8)
alpha_expr = sympy.Rational(numerator, denominator) * epsilon
alpha_decimal = float(alpha_expr)

# Print the full equation with numbers
print(f"α = ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
print(f"α = {numerator} / ({2} * {n + 8}) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")
print(f"α = {alpha_expr}")
print(f"α ≈ {alpha_decimal:.4f}")

final_answer = alpha_expr