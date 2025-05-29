import math

# Coefficients
a = -7
b = 97
c = -31

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the two solutions using the quadratic formula
if discriminant >= 0:
    sqrt_discriminant = math.sqrt(discriminant)
    n1 = (-b + sqrt_discriminant) / (2*a)
    n2 = (-b - sqrt_discriminant) / (2*a)
    # Round the solutions to 4 decimal places
    n1 = round(n1, 4)
    n2 = round(n2, 4)
    # Print the solutions
    print(f"{n1}, {n2}")
else:
    # No real solutions
    print("")