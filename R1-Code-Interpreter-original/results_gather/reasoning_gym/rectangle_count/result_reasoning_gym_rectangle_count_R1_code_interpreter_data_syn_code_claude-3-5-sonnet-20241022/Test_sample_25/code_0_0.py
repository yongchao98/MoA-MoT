def count_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [line for line in grid.split('\n') if line]
    
    height = len(grid)
    if height == 0:
        return 0
    width = len(grid[0])
    count = 0
    
    # Function to check if a position is part of a rectangle border
    def is_border(i, j):
        return 0 <= i < height and 0 <= j < width and (grid[i][j] == '#' or grid[i][j] == '█')
    
    # For each possible top-left corner
    for i in range(height):
        for j in range(width):
            if is_border(i, j):
                # Look for right edge
                right = j
                while right < width and is_border(i, right):
                    right += 1
                right -= 1
                
                if right > j:  # Found horizontal line
                    # Look for bottom edge
                    bottom = i
                    while bottom < height and is_border(bottom, j):
                        bottom += 1
                    bottom -= 1
                    
                    if bottom > i:  # Found vertical line
                        # Verify it's a complete rectangle
                        is_rectangle = True
                        # Check top and bottom edges
                        for x in range(j, right + 1):
                            if not (is_border(i, x) and is_border(bottom, x)):
                                is_rectangle = False
                                break
                        # Check left and right edges
                        for y in range(i, bottom + 1):
                            if not (is_border(y, j) and is_border(y, right)):
                                is_rectangle = False
                                break
                        
                        if is_rectangle:
                            # Skip the rest of this rectangle to avoid counting it multiple times
                            count += 1
                            j = right
    
    print(count)

# Example grid
grid = """






































































                                               #####################            
                                               #                   #            
                                               #####################            """

count_rectangles(grid)