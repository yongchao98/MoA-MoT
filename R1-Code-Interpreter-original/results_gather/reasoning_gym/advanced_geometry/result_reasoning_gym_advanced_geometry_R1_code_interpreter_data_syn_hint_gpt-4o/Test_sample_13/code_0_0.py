import math

# Coordinates of the triangle
A = (3, 5)
B = (-2, -8)
C = (-1, 0)

# Calculate the side lengths
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
CA = math.sqrt((A[0] - C[0])**2 + (A[1] - C[1])**2)

# Calculate the semi-perimeter
s = (AB + BC + CA) / 2

# Calculate the area using the coordinates
area = 0.5 * abs(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]))

# Calculate the radius of the incircle
r = area / s

# Print the radius rounded to 3 decimal places
print(round(r, 3))