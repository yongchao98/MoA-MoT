import sympy

# Define the variables
x, a, b, c, d = sympy.symbols('x a b c d')

# Based on the analysis, the numerator has roots at -b, b, and d.
# The simplest form is (x-d)*(x-b)*(x+b) = (x-d)*(x**2 - b**2).
# For aesthetic reasons, this can be written as (b**2 - x**2)*(x-d) by
# absorbing a negative sign, which helps match the asymptote behavior.
numerator = (b**2 - x**2) * (x - d)

# The denominator has vertical asymptotes at x=a and x=c.
# The simplest form is (x-a)*(x-c).
denominator = (x - a) * (x - c)

# The function f(x)
f_x = numerator / denominator

# The equation is f(x) = N(x)/D(x). Let's print the components.
# We'll expand the polynomials to show the full form.
expanded_numerator = sympy.expand(numerator)
expanded_denominator = sympy.expand(denominator)

# To present the final equation clearly, we print it term by term.
# The overall function is f(x) = [(b^2-x^2)(x-d)] / [(x-a)(x-c)]
# We can represent b^2-x^2 as (b-x)(b+x)

print("Based on the analysis, the equation for f(x) is of the form:")
print("f(x) = (N(x)) / (D(x))")
print("\nWhere the numerator N(x) is:")
# N(x) = (b^2-x^2)(x-d) = (b-x)(b+x)(x-d)
print(f"N(x) = (b - x) * (b + x) * (x - d)")

print("\nAnd the denominator D(x) is:")
# D(x) = (x-a)(x-c)
print(f"D(x) = (x - a) * (x - c)")

final_equation_str = f"f(x) = ((b**2 - x**2)*(x - d)) / ((x - a)*(x - c))"
print(f"\nSo, the final simplified equation is:\n{final_equation_str}")

# We can also show the expanded form to highlight the polynomial orders.
print("\nExpanded form:")
print(f"f(x) = ({expanded_numerator}) / ({expanded_denominator})")

# Let's extract coefficients to explicitly show each part of the equation
# Numerator: -x**3 + d*x**2 + b**2*x - b**2*d
num_coeffs = sympy.Poly(expanded_numerator, x).all_coeffs()
# Denominator: x**2 - a*x - c*x + a*c
den_coeffs = sympy.Poly(expanded_denominator, x).all_coeffs()

# This is a bit complex for the final output, let's stick to the factored form.
# The final request is to output each number in the final equation.
# Since we have variables, we print the expression structure.
print("\nFinal Equation Structure:")
print("f(x) = ( (b**2) - (x**2) ) * (x - d) / ( (x - a) * (x - c) )")
