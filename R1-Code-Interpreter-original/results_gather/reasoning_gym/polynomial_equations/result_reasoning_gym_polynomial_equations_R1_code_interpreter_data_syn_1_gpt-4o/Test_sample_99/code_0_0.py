import math

# Coefficients
a = -3
b = -14
c = -16

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Check the discriminant and calculate the roots
if discriminant > 0:
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    result = f"{root1:.4f}, {root2:.4f}"
elif discriminant == 0:
    root = -b / (2*a)
    result = f"{root:.4f}"
else:
    result = ""

print(result)