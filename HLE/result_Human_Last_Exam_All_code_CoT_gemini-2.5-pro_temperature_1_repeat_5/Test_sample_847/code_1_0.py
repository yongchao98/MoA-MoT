import numpy as np
import math

# The problem reduces to finding the cost-per-bit 'c'.
# This constant is derived from the real root of a polynomial equation.
# The equation is x^3 - x^2 - 1 = 0.
# Let's find the real root of this equation.
coeffs = [1, -1, 0, -1]
roots = np.roots(coeffs)

# The roots are complex and one is real. We need the real root.
real_root = None
for r in roots:
    if np.isreal(r):
        real_root = np.real(r)
        break

print("The final cost is determined by the real root of the equation x^3 - x^2 - 1 = 0.")
print(f"The real root (x0) is: {real_root:.6f}")

# The optimal cost per bit 'c' is given by the formula c = 1 / log2(x0).
# This can also be calculated as c = ln(2) / ln(x0).
c = math.log(2) / math.log(real_root)

print("\nThe minimal cost per bit of information is calculated as c = ln(2) / ln(x0).")
print(f"The constant factor for the minimal number of coins is: {c:.6f}")

# The total cost is c * n*log2(n) for large n. The question asks for this constant c.
# We will provide the answer rounded to 3 decimal places.
final_answer = round(c, 3)
print(f"\nThe final answer, rounded to 3 decimal places, is: {final_answer}")