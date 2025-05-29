# Number of sides
n_sides = 28

# Length of each side
side_length = 56 / n_sides

# Since each side is 2 units, the polygon can be visualized as a series of 2x2 squares
# The number of such squares is half the number of sides (since each square has 4 sides)
n_squares = n_sides // 2

# Area of each 2x2 square
area_per_square = side_length ** 2

# Total area
total_area = n_squares * area_per_square

print(total_area)