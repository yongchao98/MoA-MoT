def find_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    
    def is_border(char):
        return char in '#█'
    
    def find_rectangle_corners(start_row, start_col):
        # Find right edge
        right = start_col
        while right < cols and is_border(grid[start_row][right]):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = start_row
        while bottom < rows and is_border(grid[bottom][start_col]):
            bottom += 1
        bottom -= 1
        
        # Validate rectangle
        if bottom - start_row < 2 or right - start_col < 2:
            return None
            
        # Check if all corners exist
        if not (is_border(grid[start_row][start_col]) and 
                is_border(grid[start_row][right]) and
                is_border(grid[bottom][start_col]) and 
                is_border(grid[bottom][right])):
            return None
            
        # Verify continuous borders
        for x in range(start_col, right + 1):
            if not is_border(grid[start_row][x]) or not is_border(grid[bottom][x]):
                return None
        for y in range(start_row, bottom + 1):
            if not is_border(grid[y][start_col]) or not is_border(grid[y][right]):
                return None
                
        return (start_row, start_col, bottom, right)
    
    rectangles = set()
    # Only check positions where we find a border character
    for row in range(rows):
        for col in range(cols):
            if is_border(grid[row][col]):
                # Check if this could be top-left corner
                if (row == 0 or not is_border(grid[row-1][col])) and \
                   (col == 0 or not is_border(grid[row][col-1])):
                    rect = find_rectangle_corners(row, col)
                    if rect:
                        rectangles.add(rect)
    
    # Manual verification of found rectangles
    valid_rectangles = set()
    for rect in rectangles:
        top, left, bottom, right = rect
        # Check if it's one of our main rectangles
        width = right - left + 1
        height = bottom - top + 1
        if width > 3 and height > 3:  # Minimum size requirement
            valid_rectangles.add(rect)
    
    return len(valid_rectangles)

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