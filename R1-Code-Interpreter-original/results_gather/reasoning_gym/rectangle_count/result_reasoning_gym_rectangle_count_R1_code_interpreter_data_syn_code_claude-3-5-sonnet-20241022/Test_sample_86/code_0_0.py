def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = 0

    # Helper function to check if a rectangle exists between two corners
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        return True

    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] not in '#█':
                continue
            # Found top-left corner, now look for bottom-right corner
            for y2 in range(y1, height):
                for x2 in range(x1, width):
                    if grid[y2][x2] not in '#█':
                        continue
                    # Check if this forms a valid rectangle
                    if is_valid_rectangle(y1, x1, y2, x2):
                        rectangles += 1

    return rectangles

# Create the grid from the input
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "  #########################################                                     ",
    "  #                                       #                                     ",
    "  #                                       #                                     ",
    "  #########################################                                     ",
    "                                       #########                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                                ",
    "                                       #       #                             ###",
    "                                       #########                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             # #",
    "                                                                             ###",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))