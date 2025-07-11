import sympy

# Define symbols for the radii and pressure
a, b, P_i = sympy.symbols('a b P_i')

# The ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2)
ratio_formula = (b**2 + a**2) / (b**2 - a**2)

# We are given that the outer radius is twice the inner radius (b = 2a)
# Substitute b = 2*a into the ratio formula
ratio_calculated = ratio_formula.subs(b, 2*a)

# Simplify the expression
final_ratio = sympy.simplify(ratio_calculated)

# Let's show the calculation steps with numbers, assuming a=1, so b=2
a_val = 1
b_val = 2

numerator_val = b_val**2 + a_val**2
denominator_val = b_val**2 - a_val**2

print("Step 1: The formula for the ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2).")
print(f"Step 2: Given the outer radius 'b' is twice the inner radius 'a', we set b = 2a.")
print(f"Step 3: To calculate the ratio, let's assume an inner radius a = {a_val}. Then the outer radius b = {b_val}.")
print(f"Step 4: Substitute these values into the formula.")
print(f"   - Numerator = b^2 + a^2 = {b_val}^2 + {a_val}^2 = {b_val**2} + {a_val**2} = {numerator_val}")
print(f"   - Denominator = b^2 - a^2 = {b_val}^2 - {a_val}^2 = {b_val**2} - {a_val**2} = {denominator_val}")
print(f"Step 5: The final ratio is Numerator / Denominator.")
print(f"Final Equation: Ratio = {numerator_val} / {denominator_val}")
print(f"The calculated ratio is: {final_ratio}")
print(f"As a decimal, the ratio is approximately: {float(final_ratio):.3f}")