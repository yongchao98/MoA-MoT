def count_rectangles(grid):
    # Initialize a set to store unique rectangles
    rectangles = set()
    
    # Get the dimensions of the grid
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    
    # Function to find the bottom-right corner of a rectangle
    def find_bottom_right(r, c):
        max_r, max_c = r, c
        # Extend downwards
        while max_r + 1 < rows and grid[max_r + 1][c] in '#█':
            max_r += 1
        # Extend rightwards
        while max_c + 1 < cols and grid[r][max_c + 1] in '#█':
            max_c += 1
        return max_r, max_c
    
    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in '#█':
                # Find the bottom-right corner of the rectangle
                br_r, br_c = find_bottom_right(r, c)
                # Add the rectangle to the set
                rectangles.add((r, c, br_r, br_c))
                # Mark the rectangle as visited
                for i in range(r, br_r + 1):
                    for j in range(c, br_c + 1):
                        grid[i] = grid[i][:j] + ' ' + grid[i][j+1:]
    
    # Return the number of unique rectangles
    return len(rectangles)

# Define the grid
grid = [
    "                                                                    ####███#    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    #   # ##    ",
    "                                                                    ####█#█#    ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        # #     ",
    "                                                                        ###     ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "      ######################                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                                    ",
    "      #                    #                                      #######       ",
    "      #                    #                                      #     #       ",
    "      ######################                                      #######       ",
    "                                                                                ",
    "                                                                                ",
]

# Count the rectangles
total_rectangles = count_rectangles(grid)

# Print the result
print(total_rectangles)