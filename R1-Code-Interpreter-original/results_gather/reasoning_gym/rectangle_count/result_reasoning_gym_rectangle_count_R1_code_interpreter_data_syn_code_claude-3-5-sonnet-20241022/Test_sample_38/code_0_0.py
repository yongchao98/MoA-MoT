def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Helper function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Helper function to validate rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if is_border(grid[y1][x1]):
                # Found potential top-left corner
                for y2 in range(y1, height):
                    for x2 in range(x1, width):
                        if is_border(grid[y2][x2]):
                            # Check if this forms a valid rectangle
                            if is_valid_rectangle(y1, x1, y2, x2):
                                rectangles.append((y1, x1, y2, x2))
    
    # Remove duplicates and count unique rectangles
    unique_rectangles = set(rectangles)
    return len(unique_rectangles)

# Convert the grid into a list of strings
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
    "                                         #########                              ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))