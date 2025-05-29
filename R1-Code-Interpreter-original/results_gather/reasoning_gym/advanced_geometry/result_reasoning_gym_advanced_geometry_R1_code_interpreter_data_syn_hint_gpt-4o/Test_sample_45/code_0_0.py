import math

# Coordinates of the vertices
A = (0, -4)
B = (-3, -4)
C = (-7, -1)

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate the side lengths
AB = distance(A, B)
BC = distance(B, C)
CA = distance(C, A)

# Calculate the semi-perimeter
s = (AB + BC + CA) / 2

# Calculate the area using the coordinates
area = 0.5 * abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]))

# Calculate the incircle radius
r = area / s

# Print the incircle radius rounded to 3 decimal places
print(round(r, 3))