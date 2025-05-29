import math

# Coordinates of the vertices
A = (8, 3)
B = (4, -9)
C = (6, 5)

# Function to calculate distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate side lengths
AB = distance(A, B)
BC = distance(B, C)
CA = distance(C, A)

# Calculate semi-perimeter
s = (AB + BC + CA) / 2

# Calculate area using Heron's formula
area = math.sqrt(s * (s - AB) * (s - BC) * (s - CA))

# Calculate incircle radius
r = area / s

# Print the incircle radius rounded to 3 decimal places
print(round(r, 3))