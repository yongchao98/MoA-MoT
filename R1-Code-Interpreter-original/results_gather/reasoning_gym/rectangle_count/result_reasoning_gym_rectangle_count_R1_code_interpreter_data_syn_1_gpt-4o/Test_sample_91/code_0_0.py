# Define the grid
grid = [
    "                                  ########                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  #      #                                      ",
    "                                  ########                                      "
]

# Function to count rectangles
def count_rectangles(grid):
    # Initialize variables
    rectangle_count = 0
    visited = set()
    
    # Function to check if a point is within bounds and is a '#'
    def is_rectangle_char(x, y):
        return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] == '#'
    
    # Function to explore a rectangle starting from the top-left corner
    def explore_rectangle(x, y):
        # Find the bottom-right corner
        max_x, max_y = x, y
        while is_rectangle_char(max_x + 1, y):
            max_x += 1
        while is_rectangle_char(x, max_y + 1):
            max_y += 1
        
        # Mark all '#' in the rectangle as visited
        for i in range(x, max_x + 1):
            for j in range(y, max_y + 1):
                visited.add((i, j))
    
    # Iterate over the grid
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and grid[i][j] == '#':
                # Found a new rectangle
                rectangle_count += 1
                explore_rectangle(i, j)
    
    return rectangle_count

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)