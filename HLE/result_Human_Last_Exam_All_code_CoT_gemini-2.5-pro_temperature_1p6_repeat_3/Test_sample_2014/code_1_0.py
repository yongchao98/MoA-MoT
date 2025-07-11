import sympy

# Define symbolic variables
c, s, h = sympy.symbols('c s h')

# Given parameters
params = {
    s: c/2,
    h: c/2
}

# Step 1: Define the constants k_h and k_s
k_h_expr = c / (4 * h)
k_s_expr = c / (2 * s) - (c * s) / (2 * (s**2 + 4 * h**2))

# Step 2: Substitute the given parameters to find the numerical values
k_h_val = k_h_expr.subs(params)
k_s_val = k_s_expr.subs(params)

# Step 3: Define the formula for the lift ratio
# This formula is derived from the system of linear equations for Γ1 and Γ2
numerator_expr = 1 - k_h_val - k_s_val
denominator_expr = 1 - k_h_val + k_s_val
lift_ratio = numerator_expr / denominator_expr

# Step 4: Print the calculation steps
print("The lift ratio L1/L2 is calculated using the formula: (1 - k_h - k_s) / (1 - k_h + k_s)\n")
print(f"Given s = c/2 and h = c/2:")
print(f"k_h = c / (4 * h) = {c} / (4 * {params[h]}) = {k_h_val}")
print(f"k_s = c/(2*s) - c*s/(2*(s^2 + 4*h^2)) = {k_s_val}")
print("\nSubstituting these values into the lift ratio formula:")
print(f"L1/L2 = (1 - {k_h_val} - {k_s_val}) / (1 - {k_h_val} + {k_s_val})")
print(f"L1/L2 = ({numerator_expr}) / ({denominator_expr})")
print(f"L1/L2 = {float(lift_ratio):.4f}")
print(f"As a fraction, L1/L2 = {sympy.fraction(lift_ratio)[0]}/{sympy.fraction(lift_ratio)[1]}")

# Final answer in the required format
final_answer = float(lift_ratio)
# final_answer_str = f"<<<{final_answer:.4f}>>>"
# print(final_answer_str) # Let's output the fraction for precision in the thinking.
