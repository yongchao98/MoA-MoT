import sympy

# Define x as a symbol
x = sympy.Symbol('x')

# The function f(x) for the blue curve, derived from its features (asymptotes and extrema),
# is f(x) = x + 2 + 4/(x-2).
f_x = x + 2 + 4 / (x - 2)

# Calculate the second derivative, f''(x).
f_double_prime_x = sympy.diff(f_x, x, 2)

# Construct the target function y = -0.5 * f''(3x-2) + 1.
# Note: sympy.S(0.5) ensures it's treated as a Sympy object for precision.
y_final = -sympy.S(0.5) * f_double_prime_x.subs(x, 3*x - 2) + 1

# Simplify the final expression for y.
y_final_simplified = sympy.simplify(y_final)

# Print the final derived equation in a readable format.
final_eq_str = str(y_final_simplified).replace('**', '^')
print(f"The simplified form of the target function is: y = {final_eq_str}")
print("\nThis function has a horizontal asymptote at y=1 and a vertical asymptote at x=4/3, which matches the black curve.")

# As requested, output each number in the final equation.
print("\nThe numbers defining the original transformation y = -0.5 * f''(3*x - 2) + 1 are:")
print(f"Vertical scaling and reflection factor: -0.5")
print(f"Horizontal compression factor: 3")
print(f"Term for horizontal shift: -2")
print(f"Vertical shift: 1")