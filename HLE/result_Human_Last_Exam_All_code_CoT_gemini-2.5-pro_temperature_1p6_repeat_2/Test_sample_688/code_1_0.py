import math

# This script calculates the prefactor c_n for the fully f-connected
# Ree-Hoover diagram's contribution to the n-th virial coefficient B_n.
# The derived formula for the prefactor is c_n = -(n-1) / n!

# We will demonstrate the calculation for n=5.
n = 5

# The numerator of the fraction is -(n-1).
numerator = -(n - 1)

# The denominator is the factorial of n, n!.
denominator = math.factorial(n)

# To present the result in its simplest form, we can simplify the fraction
# by dividing the numerator and denominator by their greatest common divisor (GCD).
common_divisor = math.gcd(abs(numerator), denominator)

simplified_numerator = numerator // common_divisor
simplified_denominator = denominator // common_divisor

# Print the step-by-step calculation of the prefactor.
print(f"The general formula for the prefactor is: c_n = -(n-1) / n!")
print(f"For n = {n}, the calculation is:")
print(f"c_{n} = -({n} - 1) / {n}!")
print(f"c_{n} = {numerator} / {denominator}")

# If the fraction was simplified, print the simplified result.
if common_divisor > 1:
    print(f"c_{n} = {simplified_numerator} / {simplified_denominator}")
