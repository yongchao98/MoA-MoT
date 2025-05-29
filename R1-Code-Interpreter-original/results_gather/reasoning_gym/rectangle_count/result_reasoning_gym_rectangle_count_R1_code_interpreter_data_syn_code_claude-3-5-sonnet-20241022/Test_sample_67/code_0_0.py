def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(i, j):
        # Check if position is a corner (# or █)
        return grid[i][j] in '#█'
    
    def trace_rectangle(top, left):
        # Find bottom right corner of rectangle starting from top left
        right = left
        bottom = top
        
        # Find right edge
        while right < width and is_corner(top, right):
            right += 1
        right -= 1
        
        # Find bottom edge
        while bottom < height and is_corner(bottom, left):
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        for i in range(top, bottom + 1):
            for j in range(left, right + 1):
                if i in (top, bottom) or j in (left, right):
                    if not is_corner(i, j):
                        return None
        return (top, left, bottom, right)
    
    rectangles = set()
    # Scan for top-left corners
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):
                # Check if it's potentially a top-left corner
                if (i == 0 or not is_corner(i-1, j)) and (j == 0 or not is_corner(i, j-1)):
                    rect = trace_rectangle(i, j)
                    if rect:
                        rectangles.add(rect)
    
    # Count overlapping rectangles
    overlaps = 0
    for i in range(height):
        for j in range(width):
            if grid[i][j] == '█':
                count = sum(1 for (top, left, bottom, right) in rectangles 
                          if top <= i <= bottom and left <= j <= right)
                if count == 2:  # If two rectangles overlap here
                    overlaps += 1
    
    # Each overlap point represents two rectangles sharing that point
    return len(rectangles)

# Test grid (your ASCII grid as a string)
grid = """
                                                                                
                                                                                
[... rest of your grid ...]
                                                                                
"""

print(find_rectangles(grid))