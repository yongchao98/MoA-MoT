def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_edge(char):
        return char in '#█'
    
    def find_complete_rectangle(top, left):
        if not is_valid_edge(grid[top][left]):
            return None
            
        # Find right edge
        right = left + 1
        while right < width and is_valid_edge(grid[top][right]):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = top + 1
        while bottom < height and is_valid_edge(grid[bottom][left]):
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        # Check if right edge is complete
        for y in range(top, bottom + 1):
            if not is_valid_edge(grid[y][right]):
                return None
                
        # Check if bottom edge is complete
        for x in range(left, right + 1):
            if not is_valid_edge(grid[bottom][x]):
                return None
                
        return (top, left, bottom, right)
    
    # Find all base rectangles
    rectangles = set()
    for i in range(height):
        for j in range(width):
            rect = find_complete_rectangle(i, j)
            if rect:
                rectangles.add(rect)
    
    # Count base rectangles
    base_count = len(rectangles)
    
    # Count overlapping areas (each █ represents one additional rectangle)
    overlap_count = sum(1 for i in range(height) 
                       for j in range(width) 
                       if grid[i][j] == '█')
    
    return base_count + overlap_count

# Test with the grid
grid = [
    "   ################                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    # ... [rest of the grid remains the same]
    "        ##########################################################              "
]

print(find_rectangles(grid))