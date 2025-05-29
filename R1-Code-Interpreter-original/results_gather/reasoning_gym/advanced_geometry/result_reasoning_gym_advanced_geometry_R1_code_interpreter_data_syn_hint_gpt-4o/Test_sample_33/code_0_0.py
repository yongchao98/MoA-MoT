import math

# Coordinates of the vertices
A = (-1, -1)
B = (6, -3)
C = (4, 4)

# Calculate the side lengths
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
CA = math.sqrt((A[0] - C[0])**2 + (A[1] - C[1])**2)

# Calculate the semi-perimeter
s = (AB + BC + CA) / 2

# Calculate the area using Heron's formula
area = math.sqrt(s * (s - AB) * (s - BC) * (s - CA))

# Calculate the incircle radius
r = area / s

# Print the incircle radius rounded to 3 decimal places
print(round(r, 3))