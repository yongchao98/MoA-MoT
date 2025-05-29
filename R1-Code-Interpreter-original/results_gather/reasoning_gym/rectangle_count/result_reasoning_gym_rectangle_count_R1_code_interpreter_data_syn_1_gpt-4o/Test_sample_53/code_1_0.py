# Define the grid as a list of strings
grid = [
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    "             ############################################                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             #                                          #                       ",
    "             ############################################                       ",
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    "                                                ####################            ",
    "                                                #                  #            ",
    "                                                #                  #            ",
    "                                                #                  #            ",
    "                                                #                  #            ",
    "                                                ####################       #####",
    "                                                           ##############  #   #",
    "                                                           #            #  #   #",
    "                                                           #            #  #   #",
    "                                                           ##############  #   #",
    "                                                                           #####",
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
    " " * 80,
]

# Function to count rectangles
def count_rectangles(grid):
    rectangles = 0
    rows = len(grid)
    cols = len(grid[0])

    # Create a visited matrix to track counted rectangles
    visited = [[False] * cols for _ in range(rows)]

    # Function to check if a rectangle is closed
    def is_closed_rectangle(start_row, end_row, start_col, end_col):
        # Check top and bottom boundaries
        if not (grid[start_row][start_col:end_col+1].strip('#') == '' and
                grid[end_row][start_col:end_col+1].strip('#') == ''):
            return False
        # Check left and right boundaries
        for r in range(start_row, end_row + 1):
            if grid[r][start_col] != '#' or grid[r][end_col] != '#':
                return False
        return True

    # Iterate over each row to find rectangles
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#' and not visited[r][c]:
                # Find the end of the horizontal line
                end_col = c
                while end_col < cols and grid[r][end_col] == '#':
                    end_col += 1
                end_col -= 1

                # Find the end of the vertical line
                end_row = r
                while end_row < rows and grid[end_row][c] == '#':
                    end_row += 1
                end_row -= 1

                # Check if we have a closed rectangle
                if is_closed_rectangle(r, end_row, c, end_col):
                    rectangles += 1
                    # Mark the rectangle as visited
                    for i in range(r, end_row + 1):
                        for j in range(c, end_col + 1):
                            visited[i][j] = True

    return rectangles

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)