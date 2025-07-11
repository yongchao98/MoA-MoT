# The dimensions of the smallest integer length rectangle admitting a non-guillotine tiling
# with squares from the set S = {2x2, 3x3, 5x5, 7x7}.
width = 12
height = 13

# The number of each type of square used in the known tiling.
# Note that the 7x7 square is not used in this minimal tiling.
num_s5 = 3  # Number of 5x5 squares
num_s3 = 5  # Number of 3x3 squares
num_s2 = 9  # Number of 2x2 squares

# Calculate the area of the rectangle.
rectangle_area = width * height

# Calculate the total area of the tiles.
s5_area = num_s5 * (5**2)
s3_area = num_s3 * (3**2)
s2_area = num_s2 * (2**2)
tiles_area = s5_area + s3_area + s2_area

# The problem asks for the area of this rectangle.
# The code will print the final equation showing the calculation.
print(f"The smallest integer length rectangle has dimensions {width}x{height}.")
print(f"The tiling consists of:")
print(f"- {num_s5} squares of size 5x5, with a total area of {num_s5} * 25 = {s5_area}")
print(f"- {num_s3} squares of size 3x3, with a total area of {num_s3} * 9 = {s3_area}")
print(f"- {num_s2} squares of size 2x2, with a total area of {num_s2} * 4 = {s2_area}")
print("\nFinal Equation:")
print(f"{s5_area} + {s3_area} + {s2_area} = {tiles_area}")
print(f"The area of the rectangle is {width} * {height} = {rectangle_area}")
