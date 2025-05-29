import math

# Coordinates of A
a = (0, 0)

# Coordinates of E
angle_e = -math.pi/4 + math.atan(1/math.sqrt(2)) + math.atan(1/math.sqrt(3))
e = (2 * math.cos(angle_e), 2 * math.sin(angle_e))

# Calculate the distance EA
ea = math.sqrt((e[0] - a[0])**2 + (e[1] - a[1])**2)

# Perimeter of the pentagon
perimeter = 1 + 1 + 1 + 1 + ea

print(perimeter)