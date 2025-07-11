import numpy as np

# The problem of finding the minimal cost per bit boils down to solving
# the cubic equation x^3 + c1*x + c0 = 0, where c1=1 and c0=-1.
# The coefficients of the polynomial x^3 + 0*x^2 + 1*x - 1 are [1, 0, 1, -1].
c3, c2, c1, c0 = 1, 0, 1, -1
coeffs = [c3, c2, c1, c0]

# Find the roots of the polynomial equation.
roots = np.roots(coeffs)

# The equation has one real root and two complex conjugate roots. We need the real root.
real_root = None
for root in roots:
    if np.isreal(root):
        real_root = np.real(root)
        break

# The minimal cost per bit, lambda, is related to the real root x0 by:
# lambda = -1 / log2(x0).
# The numbers in this final equation are -1 and 2 (for the log base).
numerator = -1
log_base = 2
min_cost_per_bit = numerator / np.log2(real_root)

# Output the equation and the final answer
print(f"The equation to solve is: {c3}*x^3 + {c1}*x - {c3} = {c2}")
print(f"The real root of the equation is x0 = {real_root:.4f}")
print(f"The minimal cost per bit is calculated as: -1 / log_base_{log_base}(x0)")
print(f"The minimal cost per bit is: {min_cost_per_bit:.3f}")
