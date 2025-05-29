def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Function to verify if a rectangle exists at given coordinates
    def is_rectangle(top, left, bottom, right):
        # Check corners
        if grid[top][left] != '#' or grid[top][right] != '#' or \
           grid[bottom][left] != '#' or grid[bottom][right] != '#':
            return False
            
        # Check horizontal lines
        for c in range(left + 1, right):
            if grid[top][c] != '#' or grid[bottom][c] != '#':
                return False
                
        # Check vertical lines
        for r in range(top + 1, bottom):
            if grid[r][left] != '#' or grid[r][right] != '#':
                return False
                
        return True
    
    # Convert input string to list of lists for easier processing
    grid = [list(row) for row in grid.split('\n')]
    
    # Search for rectangles
    for top in range(rows):
        for bottom in range(top + 1, rows):
            for left in range(cols):
                for right in range(left + 1, cols):
                    if is_rectangle(top, left, bottom, right):
                        count += 1
    
    return count

# Test with the provided grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                    ############                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    #          #                                
                                    ############                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))