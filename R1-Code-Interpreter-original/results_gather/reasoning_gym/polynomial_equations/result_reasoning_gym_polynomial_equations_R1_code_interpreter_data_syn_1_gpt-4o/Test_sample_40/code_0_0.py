import math

# Coefficients
a = 40
b = 107
c = 88

# Discriminant
discriminant = b**2 - 4*a*c

# Check if the discriminant is non-negative
if discriminant >= 0:
    # Calculate the two roots
    root1 = (-b + math.sqrt(discriminant)) / (2*a)
    root2 = (-b - math.sqrt(discriminant)) / (2*a)
    # Round the roots to 4 decimal places
    root1 = round(root1, 4)
    root2 = round(root2, 4)
    # Print the roots
    print(f"0, {root1}, {root2}")
else:
    # If the discriminant is negative, there are no real roots
    print("0")