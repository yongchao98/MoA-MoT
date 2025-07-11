import numpy as np

# The equation for alpha_0 is alpha_0^6 + 2*alpha_0^4 - 1 = 0.
# Let's define the coefficients of this polynomial in alpha_0.
# The powers are 6, 4, 0.
# P(x) = 1*x^6 + 0*x^5 + 2*x^4 + 0*x^3 + 0*x^2 + 0*x - 1
coeffs = [1, 0, 2, 0, 0, 0, -1]

# Find the roots of the polynomial
roots = np.roots(coeffs)

# We are looking for the largest value, alpha_0, which must be a positive real number.
# We filter the roots to find the one that is real and positive.
# There should be only one such root.
alpha_0 = None
for r in roots:
    # Check if the root is real (imaginary part is close to zero)
    # and if it's positive.
    if np.isreal(r) and r > 0:
        alpha_0 = r.real
        break

# The final equation for alpha_0 is: a^6 + 2a^4 - 1 = 0
# The numbers in the final equation are the coefficients.
c6 = 1
c4 = 2
c0 = -1

# Print the final equation with its coefficients
print("The equation for alpha_0 is:")
print(f"{c6} * alpha_0^6 + {c4} * alpha_0^4 + {c0} = 0")

# Print the calculated value of alpha_0
print("\nThe value of alpha_0 is:")
print(alpha_0)
