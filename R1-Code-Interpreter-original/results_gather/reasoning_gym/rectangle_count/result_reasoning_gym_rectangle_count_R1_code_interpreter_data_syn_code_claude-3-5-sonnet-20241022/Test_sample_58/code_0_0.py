def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Helper function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Helper function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all possible rectangles
    for top in range(rows):
        for left in range(cols):
            if is_border(grid[top][left]):
                # Try to find bottom-right corner
                for bottom in range(top + 1, rows):
                    for right in range(left + 1, cols):
                        if (is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            if validate_rectangle(top, left, bottom, right):
                                count += 1

    return count

# Create the grid
grid = [
    " " * 80,
] * 56 + [
    "                                       #########################                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       #                       #                ",
    "                                       █##################     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █                 #     #                ",
    "                                       █       #####     #     #                ",
    "                                       █       #   #     #     #                ",
    "                                       █       #   #     #     #                ",
    "                                       █#######█###█######     #                ",
    "           ############################█#######█████###########█###########     ",
    "           #                           #                       #          #     ",
    "           #                           #                       #          #     ",
    "           #                           #                       #          #     ",
    "           #                           #########################          #     ",
    "           #                                                              #     ",
    "           ################################################################     "
]

print(find_rectangles(grid))