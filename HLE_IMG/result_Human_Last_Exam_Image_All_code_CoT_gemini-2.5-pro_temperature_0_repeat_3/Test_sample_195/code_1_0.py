import sympy

# Step 1: Define the symbolic variable
x = sympy.Symbol('x')

# Step 2: Choose simple integer values for the parameters based on the analysis
# The values satisfy the ordering -b < d < a < c < b and the relation a + c = b + d.
a = 2
b = 4
c = 3
d = 1
K = -1 # The scaling factor must be negative

# Step 3: Construct the numerator and the denominator
# Numerator from roots at -b, b, d: (x^2 - b^2)(x - d)
# Denominator from asymptotes at a, c: (x - a)(x - c)
numerator = K * (x**2 - b**2) * (x - d)
denominator = (x - a) * (x - c)

# Step 4: Expand the polynomials to get the final equation form
expanded_numerator = sympy.expand(numerator)
expanded_denominator = sympy.expand(denominator)

# Step 5: Print the final equation
# The code outputs each number (coefficient) in the final equation as requested.
print("An equation for the function f(x) is:")
print(f"f(x) = ({sympy.pretty(expanded_numerator, use_unicode=False)}) / ({sympy.pretty(expanded_denominator, use_unicode=False)})")
