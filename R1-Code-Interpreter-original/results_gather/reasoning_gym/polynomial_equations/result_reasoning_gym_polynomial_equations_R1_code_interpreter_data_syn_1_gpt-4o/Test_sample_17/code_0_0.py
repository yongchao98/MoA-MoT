import math

# Coefficients
a = 97
b = 98
c = 7

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Check the discriminant and calculate the solutions
if discriminant > 0:
    # Two real solutions
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    result = f"{root1:.4f}, {root2:.4f}"
elif discriminant == 0:
    # One real solution
    root = -b / (2*a)
    result = f"{root:.4f}"
else:
    # No real solutions
    result = ""

print(result)