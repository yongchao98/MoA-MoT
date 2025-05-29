import math

# Coordinates of the points
Y = (0, 0)
X = (36, 48)
Z = (100, 0)
W = (36, 0)

# Calculate the length of YZ
YZ = math.sqrt(60**2 + 80**2)

# Calculate the length of WZ
WZ = math.sqrt((Z[0] - W[0])**2 + (Z[1] - W[1])**2)

print(WZ)