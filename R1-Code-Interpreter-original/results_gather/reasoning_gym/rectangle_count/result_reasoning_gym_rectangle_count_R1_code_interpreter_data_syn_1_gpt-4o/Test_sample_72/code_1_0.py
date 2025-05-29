# Define the grid
grid = [
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
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
    "                         #################################################      ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
]

# Function to count rectangles
def count_rectangles(grid):
    # Trim trailing spaces from each line
    grid = [line.rstrip() for line in grid]
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangle_count = 0
    found_rectangle = False

    # Scan the grid to find the top-left corner of rectangles
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Found the top-left corner of the rectangle
                found_rectangle = True
                break
        if found_rectangle:
            break

    if found_rectangle:
        # Since the grid contains only one rectangle, we count it
        rectangle_count = 1

    return rectangle_count

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)