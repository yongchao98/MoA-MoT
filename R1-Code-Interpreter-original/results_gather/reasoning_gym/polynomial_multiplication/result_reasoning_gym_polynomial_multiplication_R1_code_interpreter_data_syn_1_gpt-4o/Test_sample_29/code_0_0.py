import math

# Coefficients
coeff1 = -588
coeff2 = -903
coeff3 = -2835

# Calculate GCD
gcd = math.gcd(math.gcd(abs(coeff1), abs(coeff2)), abs(coeff3))

print(gcd)