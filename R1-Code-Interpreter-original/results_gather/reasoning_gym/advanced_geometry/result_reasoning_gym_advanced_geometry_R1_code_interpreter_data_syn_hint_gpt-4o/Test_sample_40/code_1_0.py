# Solve for intersection with y = 7
horizontal_line_y7 = Eq(y, 7)

# Solve for intersection
intersection_y7 = solve((altitude_C, horizontal_line_y7), (x, y))
print(intersection_y7)