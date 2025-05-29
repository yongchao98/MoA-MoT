# Define the grid as a list of strings
grid = [
    "                                         #########                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                     ####█#######█###################           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     #   #       #                  #           ",
    "                                     ####█#######█###################           ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #       #                              ",
    "                                         #########                              "
]

# Function to count rectangles
def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    rectangles = 0

    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            # Check for top-left corner of a rectangle
            if grid[r][c] in '#█':
                if (r == 0 or grid[r-1][c] not in '#█') and (c == 0 or grid[r][c-1] not in '#█'):
                    # Found a top-left corner, now find the bottom-right corner
                    for br in range(r, rows):
                        if grid[br][c] not in '#█':
                            break
                        for bc in range(c, cols):
                            if grid[br][bc] not in '#█':
                                break
                            # Check if this is a bottom-right corner
                            if (br == rows-1 or grid[br+1][bc] not in '#█') and (bc == cols-1 or grid[br][bc+1] not in '#█'):
                                rectangles += 1
                                break
    return rectangles

# Count the rectangles in the grid
rectangles_count = count_rectangles(grid)
print(rectangles_count)