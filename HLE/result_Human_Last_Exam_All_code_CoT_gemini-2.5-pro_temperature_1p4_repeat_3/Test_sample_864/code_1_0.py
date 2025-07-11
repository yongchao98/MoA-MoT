import math

# Based on the spectral factorization, the whitening filter W(D) is found to be
# W(D) = (D - 1/2) / (D - 1/3)

# Define the coefficients for the numerator and denominator polynomials in D
# Numerator: a*D + b
a = 1.0
b = -1/2

# Denominator: c*D + d
c = 1.0
d = -1/3

print("The whitening filter W(D) is given by the equation:")
# We use .2f for cleaner printing of the fraction
print(f"W(D) = ({a:.1f}*D + {b:.2f}) / ({c:.1f}*D + {d:.2f})")
print("\nWhich is equivalent to:")
# For a more standard representation (powers of D^{-1} or z^{-1})
# W(D) = (1 - 0.5*D^{-1}) / (1 - 0.333*D^{-1})
num_c1 = 1.0
num_c2 = -0.5
den_c1 = 1.0
den_c2 = -1/3
print(f"W(D) = ({num_c1:.1f} + ({num_c2:.2f})*D^-1) / ({den_c1:.1f} + ({den_c2:.3f})*D^-1)")