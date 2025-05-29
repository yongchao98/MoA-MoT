# Define the grid as a list of strings
grid = [
    "                         #################################################      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #                                               #      ",
    "                         #################################################      "
]

# Initialize variables to track the number of rectangles
rectangle_count = 0

# Iterate over each line in the grid
for line in grid:
    # Check if the line contains a rectangle boundary
    if '#' in line:
        # If a rectangle boundary is found, increment the rectangle count
        rectangle_count = 1
        break

# Print the total count of rectangles
print(rectangle_count)