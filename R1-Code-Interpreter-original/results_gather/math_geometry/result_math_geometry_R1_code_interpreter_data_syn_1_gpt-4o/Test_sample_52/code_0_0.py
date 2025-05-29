# Initialize variables
side_length = 2
num_sides = 28
perimeter = 56

# Calculate the number of unit squares
# Since each side is 2 units, and the shape is a rectilinear path,
# we can simulate the path and count the number of unit squares.

# The polygon can be visualized as a series of steps or layers
# Let's assume it forms a spiral or zigzag pattern

# Calculate the number of complete squares
# Each complete square has 4 sides, so we have 28/4 = 7 complete squares
num_complete_squares = num_sides // 4

# Each complete square has an area of side_length^2
area_per_square = side_length ** 2

# Total area is the number of complete squares times the area of each square
total_area = num_complete_squares * area_per_square

print(total_area)