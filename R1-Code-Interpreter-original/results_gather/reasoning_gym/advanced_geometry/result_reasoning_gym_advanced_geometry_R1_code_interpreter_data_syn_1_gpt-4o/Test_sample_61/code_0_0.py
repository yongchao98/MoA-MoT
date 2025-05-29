import math

# Coordinates of the triangle
A = (8, 1)
B = (-7, -8)
C = (1, -6)

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate side lengths
AB = distance(A, B)
BC = distance(B, C)
CA = distance(C, A)

# Calculate semi-perimeter
s = (AB + BC + CA) / 2

# Calculate area using the determinant method
area = 0.5 * abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]))

# Calculate the radius of the incircle
r = area / s

# Print the radius rounded to 3 decimal places
print(round(r, 3))