import numpy as np

# To find the optimal cost, we solve the characteristic equation x^3 - x^2 - 1 = 0.
# The variable x is related to the cost per bit 'c' by x = 2^(1/c).
# We represent the polynomial by its coefficients: 1*x^3 - 1*x^2 + 0*x - 1.
coeffs = [1, -1, 0, -1]
roots = np.roots(coeffs)

# The equation has one real root and two complex conjugate roots. We need the real root.
real_root_mask = np.isreal(roots)
phi = roots[real_root_mask].real[0]

print("The final equation to determine the optimal strategy's cost is x^3 - x^2 - 1 = 0.")
print(f"1. The real root of this equation is phi = {phi:.10f}")

# The minimal cost per bit of information, 'c', is calculated from phi.
# From x = 2^(1/c), we have log2(x) = 1/c, so c = 1 / log2(x).
cost_per_bit = 1 / np.log2(phi)

print(f"2. The minimal cost per bit is c = 1 / log2(phi) = {cost_per_bit:.10f}")

# The question asks for the answer up to 3 decimal places.
# This represents the constant factor for the total cost for large n.
print("\nThe minimal number of coins per bit of information is:")
print(f"{cost_per_bit:.3f}")
