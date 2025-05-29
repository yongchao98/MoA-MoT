import math

# Coefficients
a = -9
b = -142
c = -27

# Discriminant
discriminant = b**2 - 4*a*c

# Calculate roots
if discriminant >= 0:
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    # Round the roots to 4 decimal places
    root1 = round(root1, 4)
    root2 = round(root2, 4)
    print(f"0, {root1}, {root2}")
else:
    print("0")