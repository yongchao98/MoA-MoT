# Given coordinates
x_A, y_A = 9, -5
x_B, y_B = 1, 7
x_C, y_C = 9, 1

# Slope of AB
slope_AB = (y_B - y_A) / (x_B - x_A)

# Perpendicular slope for altitude from C
perpendicular_slope = -1 / slope_AB

# Equation of altitude from C: y - y_C = perpendicular_slope * (x - x_C)
# We know x = 9 from the first altitude
x_orthocenter = 9
y_orthocenter = perpendicular_slope * (x_orthocenter - x_C) + y_C

# Round the result to 3 decimal places
y_orthocenter = round(y_orthocenter, 3)

# Print the orthocenter coordinates
print((x_orthocenter, y_orthocenter))