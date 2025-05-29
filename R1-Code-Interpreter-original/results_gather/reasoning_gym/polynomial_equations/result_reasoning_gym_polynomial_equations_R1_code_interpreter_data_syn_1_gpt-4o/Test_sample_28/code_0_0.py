import math

# Coefficients
a = -62
b = 6
c = 55

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Check the nature of the roots based on the discriminant
if discriminant > 0:
    # Two real solutions
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    print(f"{root1:.4f}, {root2:.4f}")
elif discriminant == 0:
    # One real solution
    root = -b / (2*a)
    print(f"{root:.4f}")
else:
    # No real solutions
    print("")