def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangle_count = 0

    # Function to check if a position is a top-left corner of a rectangle
    def is_top_left(r, c):
        if grid[r][c] != '#':
            return False
        if r > 0 and grid[r-1][c] == '#':
            return False
        if c > 0 and grid[r][c-1] == '#':
            return False
        return True

    # Function to check if a position is a bottom-right corner of a rectangle
    def is_bottom_right(r, c):
        if grid[r][c] != '#':
            return False
        if r < rows - 1 and grid[r+1][c] == '#':
            return False
        if c < cols - 1 and grid[r][c+1] == '#':
            return False
        return True

    # Traverse the grid to count rectangles
    for r in range(rows):
        for c in range(cols):
            if is_top_left(r, c):
                rectangle_count += 1
            if grid[r][c] == '█':
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                              ############      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                     #########################################█####      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #                                        #   #      #      ",
    "                     #########################################█####      #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "    ###########################################################          #      ",
    "    #                                                        ##          #      ",
    "    #                                                        ##          #      ",
    "    #                                                        ##          #      ",
    "    #                                                        ##          #      ",
    "    ###########################################################          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              #          #      ",
    "                                                              ############      "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)