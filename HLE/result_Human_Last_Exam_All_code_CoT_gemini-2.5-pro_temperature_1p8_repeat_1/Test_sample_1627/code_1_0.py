import sympy

# 1. Define the relationship between the inner and outer radii.
# The problem states the outer radius 'b' is twice the inner radius 'a'.
# We can define these symbolically to show the cancellation.
a, b = sympy.symbols('a b')
p_i = sympy.symbols('p_i')

# Relationship b = 2*a
relation = {b: 2*a}

# 2. Write the formula for the ratio of maximum tangential stress (sigma_max) to internal pressure (p_i).
# The maximum stress occurs at the inner radius (r=a).
# Ratio = (b^2 + a^2) / (b^2 - a^2)
ratio_formula = (b**2 + a**2) / (b**2 - a**2)

# 3. Substitute the given relationship into the formula.
ratio_calculated = ratio_formula.subs(relation)

# 4. Simplify the expression to find the final numerical ratio.
final_ratio = sympy.simplify(ratio_calculated)
final_ratio_float = float(final_ratio)

# 5. Print the steps and the result.
print("To find the ratio of maximum tangential stress to internal pressure, we use Lamé's equation.")
print("The ratio is given by the formula: (b² + a²) / (b² - a²)")
print("Given that the outer radius 'b' is twice the inner radius 'a' (b = 2a).")
print("\nSubstituting b = 2a into the formula:")
print("Ratio = ((2a)² + a²) / ((2a)² - a²)")
print("      = (4a² + a²) / (4a² - a²)")
print("      = 5a² / 3a²")
print("The 'a²' terms cancel out, leaving:")
print(f"Ratio = {final_ratio.p}/{final_ratio.q}")
print(f"\nAs a decimal, the ratio is approximately: {final_ratio_float:.3f}")

# Final Answer selection logic
# The calculated ratio is 5/3, which corresponds to choice E and choice H (1.667).
# Typically the exact fraction is preferred.