from sympy import symbols, sqrt, Abs

# Define the coordinates symbolically
x1, y1 = -2, 1
x2, y2 = 10, -3
x3, y3 = 5, 0

# Calculate the side lengths symbolically
AB = sqrt((x2 - x1)**2 + (y2 - y1)**2)
BC = sqrt((x3 - x2)**2 + (y3 - y2)**2)
CA = sqrt((x1 - x3)**2 + (y1 - y3)**2)

# Calculate the semi-perimeter symbolically
s = (AB + BC + CA) / 2

# Calculate the area symbolically
area = Abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / 2

# Calculate the radius of the incircle symbolically
r = area / s

# Print the radius rounded to 3 decimal places
print(round(float(r), 3))