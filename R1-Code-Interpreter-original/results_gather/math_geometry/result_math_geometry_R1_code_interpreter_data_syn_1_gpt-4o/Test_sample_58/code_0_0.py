import math

# Given values
b1 = 21
c1 = 29
c2 = 87

# Calculate the missing side a1 using the Pythagorean theorem
a1 = math.sqrt(c1**2 - b1**2)

# Calculate the scale factor
k = c2 / c1

# Calculate the shortest side of the second triangle
a2 = a1 * k

print(a2)