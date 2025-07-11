import sympy

# The problem asks for a function h(x) such that if
# -sqrt(h(b(0))) < a(0) < 0, then a(t) -> 0.
# Let x = b(0).
# Based on the analysis using a Lyapunov function, the condition for convergence
# is found to be a^2 < 2 - 2x - 4x^2.
# Therefore, h(x) = 2 - 2x - 4x^2.

# We define the function h(x) as a polynomial h(x) = c0 + c1*x + c2*x^2.
# The coefficients are:
c0 = 2
c1 = -2
c2 = -4

# We can represent the function symbolically for clear output.
x = sympy.symbols('x')
h_x = c2 * x**2 + c1 * x + c0

print("The function h(x) has been determined based on Lyapunov stability analysis.")
print(f"The function is: h(x) = {h_x}")
print("\nHere are the coefficients of the polynomial h(x) = c0 + c1*x + c2*x^2:")
print(f"The constant term, c0 = {c0}")
print(f"The coefficient of x, c1 = {c1}")
print(f"The coefficient of x^2, c2 = {c2}")
