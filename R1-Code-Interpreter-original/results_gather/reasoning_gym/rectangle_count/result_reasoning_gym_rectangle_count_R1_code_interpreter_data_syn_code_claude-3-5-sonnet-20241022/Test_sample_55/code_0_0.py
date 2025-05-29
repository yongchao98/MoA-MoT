def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    rectangles = []
    
    # Function to check if a character is part of a rectangle border
    def is_border(char):
        return char in '#█'
    
    # Function to validate if coordinates form a valid rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom borders
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right borders
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all rectangles
    for top in range(rows):
        for left in range(cols):
            if is_border(grid[top][left]):
                # Try to find bottom-right corner
                for bottom in range(top, rows):
                    for right in range(left, cols):
                        if (bottom > top and right > left and 
                            is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            if is_valid_rectangle(top, left, bottom, right):
                                rectangles.append((top, left, bottom, right))
    
    # Count unique rectangles
    return len(rectangles)

# Create the grid from the ASCII art
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
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "               ############################################################     ",
    "               #                                                          #     ",
    "               ############################################################     ",
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
    "                        ###########████████████████████████#####################",
    "                        #          #                      #                    #",
    "                        #          #                      #                    #",
    "                        ###########█######################█#####################",
    "                                   #                      #                     ",
    "                                   #                      #                     ",
    "                                   #                      #                     ",
    "                                   ########################                     ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))