import math

# Coordinates of O
O_x, O_y = intersection[x], intersection[y]

# Coordinates of A and C
A = (0, 4)
C = (4, 0)

# Calculate distances
OC = math.sqrt((O_x - C[0])**2 + (O_y - C[1])**2)
OA = math.sqrt((O_x - A[0])**2 + (O_y - A[1])**2)

# Calculate the ratio OC/OA
ratio = OC / OA
ratio