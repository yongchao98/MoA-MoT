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

def count_rectangles(grid):
    # Initialize a set to store unique rectangles
    rectangles = set()
    
    # Iterate over each row in the grid
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            # Check for the top-left corner of a rectangle
            if grid[i][j] == '#':
                # Find the bottom-right corner of the rectangle
                for k in range(i + 1, len(grid)):
                    if grid[k][j] != '#' and grid[k][j] != '█':
                        break
                    for l in range(j + 1, len(grid[k])):
                        if grid[k][l] != '#' and grid[k][l] != '█':
                            break
                        # Check if this forms a valid rectangle
                        if (grid[i][l] == '#' or grid[i][l] == '█') and (grid[k][j] == '#' or grid[k][j] == '█'):
                            # Ensure the entire boundary is valid
                            if all(grid[x][j] == '#' or grid[x][j] == '█' for x in range(i, k + 1)) and \
                               all(grid[x][l] == '#' or grid[x][l] == '█' for x in range(i, k + 1)) and \
                               all(grid[i][y] == '#' or grid[i][y] == '█' for y in range(j, l + 1)) and \
                               all(grid[k][y] == '#' or grid[k][y] == '█' for y in range(j, l + 1)):
                                rectangles.add((i, j, k, l))
    
    # Count the number of unique rectangles
    return len(rectangles)

# Count the rectangles in the grid
total_rectangles = count_rectangles(grid)
print(total_rectangles)