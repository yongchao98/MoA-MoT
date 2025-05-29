import math

# Coordinates of X, Y, Z
X = (3, 3 * math.sqrt(3))
Y = (9, 3 * math.sqrt(3))
Z = (6, 3 * math.sqrt(3) - 8)

# Distance XY
XY = math.sqrt((Y[0] - X[0])**2 + (Y[1] - X[1])**2)

# Distance YZ
YZ = math.sqrt((Z[0] - Y[0])**2 + (Z[1] - Y[1])**2)

# Distance ZX
ZX = math.sqrt((X[0] - Z[0])**2 + (X[1] - Z[1])**2)

# Perimeter of triangle XYZ
perimeter_XYZ = XY + YZ + ZX

print(perimeter_XYZ)