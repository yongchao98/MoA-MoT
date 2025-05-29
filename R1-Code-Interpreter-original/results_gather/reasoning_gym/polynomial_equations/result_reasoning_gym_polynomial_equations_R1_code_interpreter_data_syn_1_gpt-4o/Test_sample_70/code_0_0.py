import math

# Coefficients
a = -50
b = 55
c = 76

# Discriminant
discriminant = b**2 - 4*a*c

# Roots calculation
if discriminant >= 0:
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    roots = [0, root1, root2]
else:
    roots = [0]

# Round the roots to 4 decimal places
rounded_roots = [round(root, 4) for root in roots]

# Print the roots
print(rounded_roots)