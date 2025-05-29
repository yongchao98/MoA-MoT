import math

# Coefficients of the quadratic equation
a = 1
b = -6
c = 8

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Calculate the two possible solutions
x1 = (-b + math.sqrt(discriminant)) / (2*a)
x2 = (-b - math.sqrt(discriminant)) / (2*a)

# Since BP < DP, we choose the smaller value
BP = min(x1, x2)

print(BP)