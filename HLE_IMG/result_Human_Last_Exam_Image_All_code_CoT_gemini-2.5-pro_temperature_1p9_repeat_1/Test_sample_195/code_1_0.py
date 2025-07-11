import sympy

# Define the variable x and the constants a, b, c, d
x, a, b, c, d = sympy.symbols('x a b c d')

# Numerator based on x-intercepts
# Zeros at x=-b, x=b -> (x-b)(x+b) = x**2 - b**2
# Zero at x=d (touching axis) -> (x-d)**2
numerator = (x**2 - b**2) * (x - d)**2

# Denominator based on vertical asymptotes
# Asymptote at x=a (approaches +inf from both sides) -> (x-a)**2
# Asymptote at x=c (sign change) -> (x-c)
denominator = (x - a)**2 * (x - c)

# The function f(x)
# The leading coefficient is 1, as the slant asymptote is y=x
f_x = numerator / denominator

# Print the final equation in a readable format
print("The equation for the function f(x) is:")
sympy.pprint(sympy.Eq(sympy.Symbol('f(x)'), f_x), use_unicode=False)
