import math

# Coordinates of the points
A = (-9, -10)
B = (-1, -8)
C = (-10, 6)

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate the side lengths
AB = distance(A, B)
BC = distance(B, C)
CA = distance(C, A)

# Calculate the semi-perimeter
s = (AB + BC + CA) / 2

# Calculate the area using Heron's formula
area = math.sqrt(s * (s - AB) * (s - BC) * (s - CA))

# Calculate the radius of the incircle
radius = area / s

# Print the radius rounded to 3 decimal places
print(round(radius, 3))