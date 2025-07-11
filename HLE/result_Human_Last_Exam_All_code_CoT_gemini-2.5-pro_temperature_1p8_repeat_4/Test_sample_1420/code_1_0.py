import sympy

# In the context of renormalization group theory, we calculate the specific heat
# scaling exponent alpha (α) for a system in d=3 using the epsilon expansion.

# 1. Define the parameters of the problem.
d = 3  # Spatial dimension
d_c = 4  # Upper critical dimension for the O(n) model
n = 1  # For the Ising universality class

# 2. Calculate the expansion parameter, epsilon (ϵ).
# The epsilon expansion is a perturbation series in powers of ϵ = d_c - d.
epsilon = d_c - d

# 3. Use the first-order epsilon expansion formula for α.
# For the O(n) model, the formula to first order in ϵ is:
# α = (4 - n) / (2 * (n + 8)) * ϵ
# We use sympy to keep the result as a fraction for precision.
num_val = 4 - n
den_val = 2 * (n + 8)

alpha_expr_num = sympy.Integer(epsilon * num_val)
alpha_expr_den = sympy.Integer(den_val)
alpha = alpha_expr_num / alpha_expr_den

# 4. Print the calculation step-by-step.
print("Calculation of the Specific Heat Exponent α")
print("-" * 45)
print(f"Given parameters:")
print(f"  - Spatial dimension (d): {d}")
print(f"  - Upper critical dimension (d_c): {d_c}")
print(f"  - Model component (n, for Ising model): {n}\n")

print("Step 1: Calculate epsilon (ϵ = d_c - d)")
print(f"  ϵ = {d_c} - {d}")
print(f"  ϵ = {epsilon}\n")

print("Step 2: Apply the first-order expansion formula for α")
print("  Formula: α = (4 - n) / (2 * (n + 8)) * ϵ")
print(f"  α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
print(f"  α = {num_val} / (2 * {n + 8})")
print(f"  α = {num_val} / {den_val}")
print(f"  α = {alpha}\n")

# As a floating-point number for context
alpha_float = float(alpha)
print(f"The final result for the specific heat exponent α is {alpha} (approximately {alpha_float:.4f}).")