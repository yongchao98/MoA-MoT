# Calculate the range of x using the triangle inequality
d = 116  # Diameter of the circle
third_part = 80  # Third part of the journey

# Triangle inequalities
x1 = d + third_part  # 116 + 80 > x
x2 = d - third_part  # 116 + x > 80
x3 = third_part - d  # x + 80 > 116

# Calculate the valid range for x
x_min = max(0, x3)
x_max = x1

# Calculate the total distance
total_distance = d + third_part + x_min

print(total_distance)