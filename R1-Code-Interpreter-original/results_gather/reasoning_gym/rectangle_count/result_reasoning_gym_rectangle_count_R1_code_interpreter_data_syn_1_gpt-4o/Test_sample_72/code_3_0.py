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

    # Scan the grid to find the top-left corner of rectangles
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Check if this is the top-left corner of a rectangle
                if (r == 0 or grid[r-1][c] != '#') and (c == 0 or grid[r][c-1] != '#'):
                    # Check if this is a closed rectangle
                    # Find the bottom-right corner
                    bottom = r
                    while bottom < rows and grid[bottom][c] == '#':
                        bottom += 1
                    right = c
                    while right < cols and grid[r][right] == '#':
                        right += 1
                    # Verify the rectangle is closed
                    if all(grid[bottom-1][x] == '#' for x in range(c, right)) and all(grid[y][right-1] == '#' for y in range(r, bottom)):
                        rectangle_count += 1

    return rectangle_count

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)