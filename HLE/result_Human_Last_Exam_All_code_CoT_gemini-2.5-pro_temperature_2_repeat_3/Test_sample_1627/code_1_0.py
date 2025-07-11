import sympy

# Define symbolic variables for the radii and pressure.
# While we can use numbers, symbols help demonstrate the formula cancellation.
r_i, r_o, P_i = sympy.symbols('r_i r_o P_i')

# 1. State the formula for the ratio of maximum tangential stress (sigma_t_max)
#    to internal pressure (P_i) in a thick-walled cylinder.
ratio_formula = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")

# 2. Apply the given condition: the outer radius is twice the inner radius (r_o = 2 * r_i).
#    Substitute this into the formula.
ratio_with_condition = ratio_formula.subs(r_o, 2*r_i)

print("\nGiven the condition that the outer radius is twice the inner radius (r_o = 2 * r_i):")
print(f"Ratio = ((2*r_i)^2 + r_i^2) / ((2*r_i)^2 - r_i^2)")

# 3. Simplify the expression to find the numerical ratio.
simplified_ratio = sympy.simplify(ratio_with_condition)
numerator, denominator = simplified_ratio.as_numer_denom()

# 4. Print the final result, showing the numerical calculation.
# To show the numbers in the final equation as requested, let's substitute a value, e.g., r_i = 1.
# This makes r_o = 2.
r_i_val = 1
r_o_val = 2
numerator_val = r_o_val**2 + r_i_val**2
denominator_val = r_o_val**2 - r_i_val**2
final_ratio_val = numerator_val / denominator_val

print(f"\nSimplifying the terms, this becomes:")
print(f"Ratio = (4*r_i^2 + r_i^2) / (4*r_i^2 - r_i^2) = (5*r_i^2) / (3*r_i^2)")
print(f"\nThe r_i^2 terms cancel out, leaving the final equation and ratio:")
print(f"Ratio = {numerator} / {denominator}")

print(f"\nIn decimal form, this is approximately {final_ratio_val:.3f}.")