import math

# Coordinates of the vertices
A = (0, -3)
B = (-5, 2)
C = (8, -1)

# Calculate the side lengths
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
CA = math.sqrt((A[0] - C[0])**2 + (A[1] - C[1])**2)

# Calculate the semi-perimeter
s = (AB + BC + CA) / 2

# Calculate the area using Heron's formula
area = math.sqrt(s * (s - AB) * (s - BC) * (s - CA))

# Calculate the inradius
inradius = area / s

# Round the inradius to 3 decimal places
inradius_rounded = round(inradius, 3)

print(inradius_rounded)