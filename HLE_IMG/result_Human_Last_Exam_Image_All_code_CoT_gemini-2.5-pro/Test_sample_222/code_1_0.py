import numpy as np

# This script calculates the imaginary part of the sum of the contour integrals.
# The calculation is based on the residue theorem as explained above.

# The total contribution from each relevant pole to the sum inside the 2*pi*i factor is calculated.
# Total contribution = W_total(zk) * Res(f, zk)

# Contribution from the pole at z = 1.5:
# W_total(1.5) = -1
# Res(f, 1.5) = 3*sqrt(pi)/4
term1 = -1 * (3 * np.sqrt(np.pi) / 4)

# Contribution from the pole at z = -1:
# W_total(-1) = 2
# Res(f, -1) = -2/5
term2 = 2 * (-2.0 / 5.0)

# Contribution from the pole at z = -3:
# W_total(-3) = 1
# Res(f, -3) = -1/9
term3 = 1 * (-1.0 / 9.0)

# The sum of these terms is the real value multiplied by 2*pi*i to get the integral.
sum_of_terms = term1 + term2 + term3

# The imaginary part of the integral is 2*pi multiplied by this sum.
imaginary_part = 2 * np.pi * sum_of_terms

print("The final equation for the imaginary part is:")
# Printing each number in the final equation as requested.
print(f"Im(Sum) = 2 * pi * ( {term1} + {term2} + {term3} )")

print("\nThe numerical result is:")
print(imaginary_part)