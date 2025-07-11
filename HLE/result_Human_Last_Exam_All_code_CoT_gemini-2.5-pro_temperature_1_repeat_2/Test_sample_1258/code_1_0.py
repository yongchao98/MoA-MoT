# This script provides the analytical expression for the fluxmetric demagnetizing factor (N_f)
# for a magnetic cylinder based on the provided definitions.

# The definitions are:
# g = length-to-diameter ratio
# k^2 = 1 / (1 + g^2 / 4)
# F(k) = complete elliptic integral of the first kind with modulus k
# E(k) = complete elliptic integral of the second kind with modulus k

# The derivation starts from the formula by R. I. Joseph (J. Appl. Phys. 37, 4639, 1966),
# which uses a modulus definition identical to the one provided:
# N_f = 1 - (2 / (g * pi * sqrt(1 + g^2/4))) * (F(k) - E(k))

# By rearranging the definition of k, we find g in terms of k:
# g = (2 * sqrt(1 - k^2)) / k

# Substituting g back into the formula for N_f and simplifying yields the final expression.

print("The analytical expression for the fluxmetric demagnetizing factor (N_f) is:")
# The final equation is printed below, showing each number and symbol.
print("N_f = 1 - (k^2 / (pi * sqrt(1 - k^2))) * (F(k) - E(k))")

print("\nWhere:")
print("k is the modulus, defined as k^2 = 1 / (1 + g^2 / 4), with g being the length-to-diameter ratio.")
print("F(k) is the complete elliptic integral of the first kind.")
print("E(k) is the complete elliptic integral of the second kind.")
print("pi is the mathematical constant pi.")
print("sqrt() denotes the square root function.")
