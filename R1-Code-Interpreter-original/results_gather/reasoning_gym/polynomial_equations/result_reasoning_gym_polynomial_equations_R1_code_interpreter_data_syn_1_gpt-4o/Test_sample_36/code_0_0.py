import math

# Coefficients
a = 115
b = 0
c = -43

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the two solutions
z1 = (-b + math.sqrt(discriminant)) / (2*a)
z2 = (-b - math.sqrt(discriminant)) / (2*a)

# Round the solutions to four decimal places
z1 = round(z1, 4)
z2 = round(z2, 4)

print(f"{z1}, {z2}")