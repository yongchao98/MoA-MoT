import math

# This script calculates the constant factor in the asymptotic formula for the minimal cost.

def solve_cubic_equation():
    """
    Numerically solves the equation p^3 + p - 1 = 0 for p in (0, 1)
    using the bisection method.
    """
    low = 0.0
    high = 1.0
    # 100 iterations are more than enough for high precision.
    for _ in range(100):
        p0 = (low + high) / 2
        if p0**3 + p0 - 1 < 0:
            low = p0
        else:
            high = p0
    return (low + high) / 2

# The asymptotic minimal cost to sort the array is C * n * ln(n) for large n.
# We need to find the constant C.
# First, we find the optimal probability p0.
p0 = solve_cubic_equation()

# The constant C is derived from p0 as C = 1 / (-ln(p0)).
# Here, math.log represents the natural logarithm (ln).
constant_C = -1 / math.log(p0)

# As per the problem's request, we print the number in the final equation,
# which is the constant C, rounded to three decimal places.
# The final equation is: Cost ≈ 2.616 * n * ln(n)
print(f"The constant C in the equation Cost ≈ C * n * ln(n) is:")
print(f"{constant_C:.3f}")
