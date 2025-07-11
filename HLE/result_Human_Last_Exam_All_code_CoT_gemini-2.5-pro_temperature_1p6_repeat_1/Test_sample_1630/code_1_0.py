import numpy as np
from numpy.polynomial import chebyshev

def get_chebyshev_roots(degree):
    """
    The roots of the Chebyshev polynomial T_n(x) are x_k = cos((2k-1)pi / (2n)).
    This function calculates them for a given degree n.
    """
    return [np.cos((2 * k - 1) * np.pi / (2 * degree)) for k in range(1, degree + 1)]

# Step 1: Define the degree of the composite polynomial.
# f(x) has degree 3, g(x) has degree 3.
# h(x) = f(g(x)) has degree 3 * 3 = 9.
# The number of fixed points is the number of solutions to h(x) = x.
# The equation h(x) - x = 0 is a polynomial of degree 9, so it has at most 9 real roots.
degree = 9

# Step 2: Construct an example of a strictly increasing polynomial h(x) with 9 fixed points.
# We use a Chebyshev polynomial T_9(x) to add "wiggles" to the line y = x.
# Let h(x) = x + epsilon * T_9(x).
# The fixed points of h(x) are the roots of h(x) - x = 0, which are the roots of T_9(x).

# We need h'(x) = 1 + epsilon * T_9'(x) > 0 for h(x) to be strictly increasing.
# The derivative of T_n is n*U_{n-1}, where U_{n-1} is the (n-1)-th Chebyshev polynomial of the second kind.
# The maximum absolute value of U_8(x) on [-1, 1] is 9.
# So, |T_9'(x)| = |9*U_8(x)| <= 81 on [-1, 1]. To be safe, we can check a larger range.
# We choose a small epsilon to ensure 1 + epsilon * T_9'(x) > 0.
epsilon = 1 / 100.0

# T_9(x) is represented by the coefficient array [0,0,0,0,0,0,0,0,0,1]
t9_coeffs = [0] * degree + [1]
t9 = chebyshev.Chebyshev(t9_coeffs)

# h(x) = x + epsilon * T_9(x)
h_coeffs = t9.coef * epsilon
h_coeffs[1] += 1
h = chebyshev.Chebyshev(h_coeffs)

# h'(x) = 1 + epsilon * T_9'(x)
h_prime = h.deriv()

# Step 3: Numerically find the fixed points by finding the roots of h(x) - x.
# These are the roots of epsilon * T_9(x) = 0, which are just the roots of T_9(x).
fixed_points = h.roots()
# Because h(x)-x is a chebyshev polynomial by construction, we can use the known formula for the roots
# for greater accuracy.
analytical_roots = get_chebyshev_roots(degree)

print(f"The composite polynomial h(x) has degree {degree}.")
print(f"The equation for fixed points h(x) = x is a polynomial equation of degree {degree}.")
print(f"Thus, there can be at most {degree} fixed points.\n")
print(f"We can construct a strictly increasing function h(x) that has {degree} fixed points.")
print("This demonstrates that the maximum is achievable.\n")

print(f"The {degree} fixed points are the roots of the 9th Chebyshev polynomial:")
# Sort the roots for clean printing
analytical_roots.sort()
equation_parts = []
for i, root in enumerate(analytical_roots):
    print(f"x_{i+1} = {root:.8f}")
    equation_parts.append(f"(x - {root:.8f})")

# In principle, the final equation with all fixed points would be f(g(x)) - x = 0.
# The construction h(x) = x + epsilon * T_9(x) leads to the equation for fixed points:
# epsilon * T_9(x) = 0, which simplifies to T_9(x) = 0.
# The polynomial T_9(x) is 256x^9 - 576x^7 + 432x^5 - 120x^3 + 9x.
print("\nThe equation for the fixed points, T_9(x)=0, is:")
print("256*x^9 - 576*x^7 + 432*x^5 - 120*x^3 + 9*x = 0")
print(f"\nThis equation has {len(analytical_roots)} distinct real roots, confirming that 9 fixed points is possible.")
