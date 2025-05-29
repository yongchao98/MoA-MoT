def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    # Find all rectangles
    rectangles = []
    height = len(grid)
    width = len(grid[0])
    
    # Function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Function to validate rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom borders
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right borders
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if is_border(grid[y1][x1]):
                # Found potential top-left corner
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if is_border(grid[y2][x2]):
                            # Check if this forms a valid rectangle
                            if is_valid_rectangle(y1, x1, y2, x2):
                                rectangles.append((y1, x1, y2, x2))
    
    # Count overlapping rectangles
    overlaps = set()
    for i, rect1 in enumerate(rectangles):
        for j, rect2 in enumerate(rectangles[i+1:], i+1):
            y1a, x1a, y2a, x2a = rect1
            y1b, x1b, y2b, x2b = rect2
            
            # Check for overlap
            if not (x2a < x1b or x1a > x2b or y2a < y1b or y1a > y2b):
                # Found overlap, mark both rectangles
                overlaps.add(i)
                overlaps.add(j)
    
    # Calculate total rectangles
    # Each '█' represents 2 overlapping rectangles
    total_rectangles = len(rectangles)
    
    print(total_rectangles)

# Create the grid string
grid = """
[Your ASCII grid here]
"""

find_rectangles(grid)