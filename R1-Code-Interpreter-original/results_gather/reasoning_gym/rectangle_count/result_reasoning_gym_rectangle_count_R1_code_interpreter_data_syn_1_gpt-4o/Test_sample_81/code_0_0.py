def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangle_count = 0

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == '#':
                # Check for the bottom-right corner of the rectangle
                for k in range(i + 1, rows):
                    if grid[k][j] != '#':
                        break
                    for l in range(j + 1, cols):
                        if grid[i][l] != '#':
                            break
                        if grid[k][l] == '#':
                            # Check if the sides are complete
                            if all(grid[m][j] == '#' for m in range(i, k + 1)) and \
                               all(grid[m][l] == '#' for m in range(i, k + 1)) and \
                               all(grid[i][n] == '#' for n in range(j, l + 1)) and \
                               all(grid[k][n] == '#' for n in range(j, l + 1)):
                                rectangle_count += 1
                                # Mark the rectangle as counted
                                for m in range(i, k + 1):
                                    for n in range(j, l + 1):
                                        grid[m][n] = '.'
    return rectangle_count

# Define the grid
grid = [
    "                                   ####                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   ####                                         "
]

# Convert grid to a list of lists
grid = [list(row) for row in grid]

# Count rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)