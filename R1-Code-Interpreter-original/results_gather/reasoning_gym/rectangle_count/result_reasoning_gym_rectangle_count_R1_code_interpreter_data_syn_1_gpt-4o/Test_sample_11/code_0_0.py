def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangle_count = 0

    # Function to check if a rectangle starts at (r, c)
    def is_rectangle_start(r, c):
        if grid[r][c] != '#':
            return False
        # Check if it's the top-left corner of a rectangle
        if (r == 0 or grid[r-1][c] == ' ') and (c == 0 or grid[r][c-1] == ' '):
            return True
        return False

    # Scan the grid
    for r in range(rows):
        for c in range(cols):
            if is_rectangle_start(r, c):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                    ############                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    ############                                "
]

# Count rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)