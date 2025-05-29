import math

# Coordinates of point E
r_E = 2
theta_E = -math.pi/4 + math.atan(1/math.sqrt(2)) + math.atan(1/math.sqrt(3))
x_E = r_E * math.cos(theta_E)
y_E = r_E * math.sin(theta_E)

# Coordinates of point A
x_A, y_A = 0, 0

# Calculate the distance EA
EA = math.sqrt((x_E - x_A)**2 + (y_E - y_A)**2)

# Calculate the perimeter
perimeter = 1 + 1 + 1 + 1 + EA

print(perimeter)