import math

# Coefficients
a = 50
b = 0
c = -89

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the two solutions
o1 = (-b + math.sqrt(discriminant)) / (2*a)
o2 = (-b - math.sqrt(discriminant)) / (2*a)

# Round the solutions to 4 decimal places
o1 = round(o1, 4)
o2 = round(o2, 4)

print(o1, o2)