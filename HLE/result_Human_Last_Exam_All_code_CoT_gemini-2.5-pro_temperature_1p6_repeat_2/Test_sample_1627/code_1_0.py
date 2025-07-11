import sympy

# Define symbols for the radii and the ratio
a, b = sympy.symbols('a b')

# The problem states that the outer radius (b) is twice the inner radius (a)
# We can represent this relationship by setting b = 2a
b_val = 2 * a

# Formula for the ratio of maximum tangential stress to internal pressure
# Ratio = (b^2 + a^2) / (b^2 - a^2)
ratio_formula = (b**2 + a**2) / (b**2 - a**2)

# Substitute b = 2a into the formula
final_ratio = ratio_formula.subs(b, b_val)

# Simplify the expression to get the numerical value
# The 'a' symbol will cancel out
simplified_ratio = sympy.simplify(final_ratio)

# To display the calculation steps, let's use numerical values.
# The ratio is independent of the actual radius, so we can set a=1 and b=2.
a_num = 1
b_num = 2

numerator = b_num**2 + a_num**2
denominator = b_num**2 - a_num**2

print(f"The formula for the ratio of maximum tangential stress to internal pressure is:")
print(f"Ratio = (b² + a²)/(b² - a²)")
print(f"Given that the outer radius is twice the inner radius (b = 2a), we can substitute values like a=1 and b=2.")
print(f"Ratio = ({b_num}² + {a_num}²)/({b_num}² - {a_num}²)")
print(f"Ratio = ({b_num**2} + {a_num**2})/({b_num**2} - {a_num**2})")
print(f"Ratio = {numerator}/{denominator}")
print(f"The final ratio is {simplified_ratio}.")
print(f"As a decimal, this is approximately {float(simplified_ratio):.3f}.")