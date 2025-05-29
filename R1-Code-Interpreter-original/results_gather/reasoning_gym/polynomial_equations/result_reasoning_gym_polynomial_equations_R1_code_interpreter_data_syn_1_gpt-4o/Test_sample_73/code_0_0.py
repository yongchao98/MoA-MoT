import math

# Coefficients
a = 68
b = -19
c = -88

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the roots if the discriminant is non-negative
if discriminant >= 0:
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    # Round the roots to 4 decimal places
    root1 = round(root1, 4)
    root2 = round(root2, 4)
    print(f"{root1}, {root2}")
else:
    # No real roots
    print("")