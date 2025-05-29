def count_rectangles(grid_str):
    # Split into lines and remove any empty lines at the start/end
    grid = [line for line in grid_str.splitlines() if line]
    height = len(grid)
    width = len(grid[0]) if grid else 0
    
    # Count distinct rectangles and overlaps
    rectangles = 0
    overlaps = 0
    
    # Helper function to validate if a point is part of rectangle border
    def is_border(y, x):
        return 0 <= y < height and 0 <= x < width and (grid[y][x] in '#█')
    
    # Helper function to verify if a point is top-left corner of a rectangle
    def is_top_left_corner(y, x):
        if not is_border(y, x):
            return False
        # Must have space above and to the left
        above_clear = y == 0 or not is_border(y-1, x)
        left_clear = x == 0 or not is_border(y, x-1)
        return above_clear and left_clear
    
    # For each potential top-left corner
    for y in range(height):
        for x in range(width):
            if is_top_left_corner(y, x):
                # Verify it's actually a rectangle by finding right and bottom borders
                cur_x = x
                cur_y = y
                
                # Find right edge
                while cur_x < width and is_border(y, cur_x):
                    cur_x += 1
                
                # Find bottom edge
                while cur_y < height and is_border(cur_y, x):
                    cur_y += 1
                
                # Verify all borders exist
                is_valid = True
                for i in range(y, cur_y):
                    if not (is_border(i, x) and is_border(i, cur_x-1)):
                        is_valid = False
                        break
                for j in range(x, cur_x):
                    if not (is_border(y, j) and is_border(cur_y-1, j)):
                        is_valid = False
                        break
                
                if is_valid:
                    rectangles += 1
            
            # Count overlaps
            if grid[y][x] == '█':
                overlaps += 1
    
    # Each overlap point represents two rectangles, so we need to add one more
    # to account for the second rectangle
    total = rectangles + overlaps
    
    print(f"Found {rectangles} base rectangles and {overlaps} overlaps")
    print(f"Total: {total}")

# [Same grid string as before]