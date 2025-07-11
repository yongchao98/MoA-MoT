import sympy

# Define the symbols
x, l, n, E0 = sympy.symbols('x l n E(0)')

# The derived expression for the kinetic energy E as a function of extension x
# E(x) = E(0) * exp(x**2 / (n**2 * l**2))
E_x = E0 * sympy.exp(x**2 / (n**2 * l**2))

# The force F is the negative derivative of the energy E with respect to x
# F = -dE/dx
F = -sympy.diff(E_x, x)

# Simplify the resulting expression
F_simplified = sympy.simplify(F)

# Extract the numbers from the simplified expression for the final print statement
# The expression is of the form: -2 * E(0) * x * exp(...) / (l**2 * n**2)
# The numbers are -2, 2, 2 in the expression and 2 in the exponential.
# For clarity in the final printout, we will construct the string manually.

force_law_str = f"F(x) = - (2 * E(0) * x / (n**2 * l**2)) * exp(x**2 / (n**2 * l**2))"
print("The derived force law for the thermally isolated polymer is:")
print(force_law_str)

# For small x, we can approximate exp(y) ≈ 1 + y. Let's show this as well.
y = x**2 / (n**2 * l**2)
F_approx = - (2 * E0 * x / (n**2 * l**2)) * (1 + y)
F_approx_simplified = sympy.simplify(F_approx)

print("\nFor small x, the force law can be approximated by a linear restoring force (Hooke's Law):")
# The leading term of the approximation
F_linear = - (2 * E0 * x) / (n**2 * l**2)
print(f"F(x) ≈ - (2 * E(0) * x) / (n**2 * l**2)")
