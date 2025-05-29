from sympy import symbols, sqrt, Abs

# Define the coordinates symbolically
x1, y1 = symbols('x1 y1')
x2, y2 = symbols('x2 y2')
x3, y3 = symbols('x3 y3')

# Assign the actual values to the symbols
coords = {x1: -2, y1: 1, x2: 10, y2: -3, x3: 5, y3: 0}

# Calculate the side lengths symbolically
AB = sqrt((x2 - x1)**2 + (y2 - y1)**2).subs(coords)
BC = sqrt((x3 - x2)**2 + (y3 - y2)**2).subs(coords)
CA = sqrt((x1 - x3)**2 + (y1 - y3)**2).subs(coords)

# Calculate the semi-perimeter symbolically
s = (AB + BC + CA) / 2

# Calculate the area symbolically
area = Abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / 2
area = area.subs(coords)

# Calculate the radius of the incircle symbolically
r = area / s

# Print the radius rounded to 3 decimal places
print(round(float(r), 3))