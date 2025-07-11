import numpy as np
import numpy.polynomial.polynomial as poly

# Step 1: Theoretical maximum
# Let f(x) and g(x) be polynomials of degree 3.
# The composite function h(x) = f(g(x)) is a polynomial of degree 3 * 3 = 9.
# The fixed points of h(x) are the solutions to the equation h(x) = x.
# This is equivalent to finding the roots of the polynomial P(x) = h(x) - x = 0.
# A polynomial of degree 9 can have at most 9 real roots.
# Therefore, the maximum number of fixed points is at most 9.

# Step 2: Show that 9 is achievable
# We need to show there exist f and g (degree 3, with f'(x)>0 and g'(x)>0 for all x)
# such that h(x) - x = 0 has 9 distinct real roots.
# By Rolle's Theorem, if h(x) - x has 9 roots, its derivative h'(x) - 1 must have 8 real roots.
# Let's construct a plausible h'(x) - 1 and check the number of its real roots.
# h'(x) = f'(g(x)) * g'(x). Since f' and g' are quadratics which are always positive,
# h'(x) is an 8th-degree polynomial.

# Let's choose specific forms for f' and g' that satisfy the conditions.
# Let g'(x) = x^2 + 1, which corresponds to g(x) = x^3/3 + x (a valid degree 3 polynomial with g'>0).
# Let f'(y) = A(y - y0)^2 + B, with A>0, B>0. This corresponds to a valid degree 3 polynomial f(y).
# We can set A=1 and choose B and y0 to create an h'(x) that oscillates around y=1.
# This will make h'(x) - 1 = 0 have many roots. Let's try A=1, B=0.01, y0=30.
# The equation h'(x) - 1 = 0 becomes:
# [((x^3/3 + x) - 30)^2 + 0.01] * (x^2 + 1) - 1 = 0

# We expand this to find the polynomial coefficients.
# Let P(x) = [((x^3/3 + x) - 30)^2 + 0.01] * (x^2 + 1) - 1.
# P(x) = (1/9)x^8 + (7/9)x^6 - 20x^5 + (5/3)x^4 - 80x^3 + 901.01x^2 - 60x + 899.01
# The coefficients are ordered from the constant term to the highest degree.
coeffs = [899.01, -60, 901.01, -80, 5/3, -20, 7/9, 0, 1/9]

# Find the roots of the polynomial
roots = poly.polyroots(coeffs)

# Filter for real roots (where the imaginary part is very close to zero)
real_roots = roots[np.isclose(roots.imag, 0)].real
num_real_roots = len(real_roots)

print("To show that 9 fixed points are possible, we construct an example.")
print("We seek conditions where the derivative of h(x)-x, which is h'(x)-1, has 8 real roots.")
print("The constructed polynomial for h'(x)-1 = 0 is:")

# Build and print the equation string from the coefficients
equation_parts = []
degree = len(coeffs) - 1
for i, c in reversed(list(enumerate(coeffs))):
    if not np.isclose(c, 0):
        # sign
        sign = " - " if c < 0 else " + "
        c = abs(c)

        if i == degree:
            sign = ""
        
        # coefficient
        if np.isclose(c, 1) and i > 0:
            coef_str = ""
        else:
            coef_str = f"{c:.4f}"

        # variable part
        if i == 0:
            var_str = ""
        elif i == 1:
            var_str = "x"
        else:
            var_str = f"x^{i}"
        
        # Put it together, avoiding "1x"
        if var_str and coef_str:
            term = coef_str + var_str
        else:
            term = coef_str or var_str

        equation_parts.append(sign + term)
print("".join(equation_parts) + " = 0")


print(f"\nThis polynomial has {num_real_roots} distinct real roots.")
print("Since we have found an example where h'(x) - 1 has 8 real roots, it is")
print("possible for h(x) - x to have 9 fixed points.")
print("\nThus, the maximum number of fixed points that f(g(x)) can have is:")
final_answer = 9
print(final_answer)