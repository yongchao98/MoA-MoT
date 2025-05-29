# Define the grid as a list of strings
grid = [
    "                                                                ############### ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                                #             # ",
    "                                                          ######█##########   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     #         #   # ",
    "                                                          #     ##########█#### ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                  ########█###############█#    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       #               ##    ",
    "                                                  #       ##################    ",
    "                                                  #                        #    ",
    "                                                  #                        #    ",
    "                                                  #                        #    ",
    "                                                  ##########################    "
]

def count_rectangles(grid):
    # Initialize variables
    rectangles = 0
    visited = set()
    
    # Function to check if a position is within bounds and is a rectangle boundary
    def is_boundary(x, y):
        return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and (grid[x][y] == '#' or grid[x][y] == '█')
    
    # Function to explore a rectangle starting from a top-left corner
    def explore_rectangle(x, y):
        # Find the bottom-right corner
        max_x, max_y = x, y
        while is_boundary(max_x + 1, y):
            max_x += 1
        while is_boundary(x, max_y + 1):
            max_y += 1
        
        # Mark all positions of this rectangle as visited
        for i in range(x, max_x + 1):
            for j in range(y, max_y + 1):
                visited.add((i, j))
    
    # Iterate over the grid to find rectangles
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and is_boundary(i, j):
                # Found a new rectangle
                rectangles += 1
                explore_rectangle(i, j)
    
    return rectangles

# Count the rectangles in the grid
total_rectangles = count_rectangles(grid)
print(total_rectangles)