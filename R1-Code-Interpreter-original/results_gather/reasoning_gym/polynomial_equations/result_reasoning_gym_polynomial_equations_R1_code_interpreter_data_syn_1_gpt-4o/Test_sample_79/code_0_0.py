import math

# Coefficients
a = 54
b = -90
c = 22

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the two roots using the quadratic formula
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