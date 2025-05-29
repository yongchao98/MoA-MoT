import math

# Given values
R = 12
h = 24

# Calculate the radius of the inscribed sphere
r = (R * h) / (math.sqrt(R**2 + h**2) + R)

# Express r in the form a*sqrt(c) - a
# We know r = a(sqrt(c) - 1), so a = r / (sqrt(c) - 1)
# Let's calculate the values of a and c
a = r / (math.sqrt(5) - 1)
c = 5

# Calculate a + c
result = a + c

print(result)