def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_corner(i, j):
        return grid[i][j] in '#█'
    
    def verify_rectangle(top, left, bottom, right):
        # Check if it forms a valid rectangle
        for i in range(top, bottom + 1):
            for j in range(left, right + 1):
                if i in (top, bottom) or j in (left, right):
                    if not grid[i][j] in '#█':
                        return False
        return True
    
    rectangles = set()
    # Find all possible rectangles
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):  # Found potential top-left corner
                # Look for bottom-right corner
                for bottom in range(i, height):
                    for right in range(j, width):
                        if is_corner(bottom, right):
                            if verify_rectangle(i, j, bottom, right):
                                rectangles.add((i, j, bottom, right))
    
    # Count overlapping rectangles
    overlaps = 0
    for i in range(height):
        for j in range(width):
            if grid[i][j] == '█':
                overlaps += 1
                
    # Each '█' character represents an overlap of 2 rectangles
    # So we need to subtract overlaps from total found rectangles
    total_rectangles = len(rectangles) - overlaps
    
    print(total_rectangles)

# Create the grid string
grid = """
[Your ASCII grid here]
"""

find_rectangles(grid)