def count_rectangles(grid):
    # Convert the grid into a list of strings
    grid = [line for line in grid.split('\n') if line]
    
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    count = 0
    
    # Function to check if a point is part of a rectangle border
    def is_border(char):
        return char in ['#', '█']
    
    # Function to validate if coordinates are within grid
    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols
    
    # Function to check if we found a complete rectangle
    def verify_rectangle(top, left, bottom, right):
        # Check horizontal borders
        for c in range(left, right + 1):
            if not is_border(grid[top][c]) or not is_border(grid[bottom][c]):
                return False
        # Check vertical borders
        for r in range(top, bottom + 1):
            if not is_border(grid[r][left]) or not is_border(grid[r][right]):
                return False
        return True
    
    # Find top-left corners and verify rectangles
    for i in range(rows):
        for j in range(cols):
            if is_border(grid[i][j]):
                # Find possible bottom-right corners
                for bottom in range(i + 1, rows):
                    for right in range(j + 1, cols):
                        if is_border(grid[bottom][right]):
                            if verify_rectangle(i, j, bottom, right):
                                count += 1
    
    return count

# Test grid
grid = """






                                                                          ###   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                                                                          ###   





































































"""

print(count_rectangles(grid))