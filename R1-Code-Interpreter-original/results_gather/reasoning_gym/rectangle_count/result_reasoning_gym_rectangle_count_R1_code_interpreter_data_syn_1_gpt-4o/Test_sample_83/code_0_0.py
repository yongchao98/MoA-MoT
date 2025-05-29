# Define the grid as a list of strings
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
    "                         ####################################################   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                            ##### #   ",
    "                         #                                            #   # #   ",
    "                         #                                            #   # #   ",
    "                         #                                            #   # #   ",
    "                         #                                            ##### #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         #                                                  #   ",
    "                         ####################################################   ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                          ######################################",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                    ######█#################                   #",
    "                                    █#####█#######         #                   #",
    "                                    █     #      #         #                   #",
    "                                    █#####█######█##########                   #",
    "                                    #     #      #                             #",
    "                                    #     #      #                             #",
    "                                    #     #      #                             #",
    "                                    #     #      #                             #",
    "                                    #     #      #                             #",
    "                                    #     #      #                             #",
    "                                    ######█#######                             #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          #                                    #",
    "                                          ######################################"
]

def count_rectangles(grid):
    # Initialize variables
    rectangles = 0
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0

    # Scan the grid to find rectangles
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Check for top-left corner of a rectangle
                if (r == 0 or grid[r-1][c] != '#') and (c == 0 or grid[r][c-1] != '#'):
                    # Find the bottom-right corner
                    br_r, br_c = r, c
                    while br_r < rows and grid[br_r][c] == '#':
                        br_r += 1
                    while br_c < cols and grid[r][br_c] == '#':
                        br_c += 1
                    # Count this rectangle
                    rectangles += 1
            elif grid[r][c] == '█':
                # Count overlapping rectangles
                rectangles += 1

    return rectangles

# Count the rectangles in the grid
rectangles_count = count_rectangles(grid)
print(rectangles_count)