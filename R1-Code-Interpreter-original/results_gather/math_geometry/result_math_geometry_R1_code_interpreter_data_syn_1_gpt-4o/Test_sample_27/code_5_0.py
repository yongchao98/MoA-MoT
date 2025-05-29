import math

# Given values
R = 12
h = 24

# Calculate the radius of the inscribed sphere
r = (R * h) / (math.sqrt(R**2 + h**2) + R)

# Express r in the form a*sqrt(c) - a
# We need to find a and c such that r = a(sqrt(c) - 1)

# Calculate the components
# Let's assume a = R / (sqrt(R^2 + h^2) + R)
# and c = (R^2 + h^2) / R^2

# Calculate the components
a = R / (math.sqrt(R**2 + h**2) + R)
c = (R**2 + h**2) / R**2

# Calculate a + c
result = a + c

print(result)