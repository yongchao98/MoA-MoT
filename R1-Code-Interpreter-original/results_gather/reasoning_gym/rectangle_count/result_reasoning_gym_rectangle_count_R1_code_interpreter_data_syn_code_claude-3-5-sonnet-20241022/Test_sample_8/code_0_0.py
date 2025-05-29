def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = 0
    
    # Helper function to check if a point is part of rectangle border
    def is_border(char):
        return char in ['#', '█']
    
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
    for top in range(height):
        for left in range(width):
            if is_border(grid[top][left]):
                # Find possible bottom and right coordinates
                for bottom in range(top, height):
                    for right in range(left, width):
                        if (is_border(grid[top][left]) and 
                            is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            
                            if validate_rectangle(top, left, bottom, right):
                                # Count overlapping sections
                                overlap_count = 0
                                for y in range(top, bottom + 1):
                                    for x in range(left, right + 1):
                                        if grid[y][x] == '█':
                                            overlap_count += 1
                                
                                # Add base rectangle
                                if overlap_count > 0:
                                    rectangles += 2  # Overlapping means 2 rectangles
                                else:
                                    rectangles += 1

    # Adjust for double counting of overlapping rectangles
    return rectangles // 2

# Parse the input grid
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

# Count rectangles
result = find_rectangles(grid)
print(f"<<<{result}>>>")