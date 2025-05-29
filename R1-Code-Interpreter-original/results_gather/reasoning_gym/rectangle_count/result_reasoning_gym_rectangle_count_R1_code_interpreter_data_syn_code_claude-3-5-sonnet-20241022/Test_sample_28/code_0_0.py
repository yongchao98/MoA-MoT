def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    def is_border(char):
        return char in ['#', '█']
    
    def verify_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    def find_overlaps():
        overlaps = set()
        for y in range(height):
            for x in range(width):
                if grid[y][x] == '█':
                    overlaps.add((y, x))
        return overlaps
    
    # Find all potential rectangles
    for top in range(height):
        for left in range(width):
            if is_border(grid[top][left]):
                # Find potential bottom and right boundaries
                for bottom in range(top, height):
                    for right in range(left, width):
                        if (is_border(grid[top][left]) and 
                            is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            if verify_rectangle(top, left, bottom, right):
                                rectangles.append((top, left, bottom, right))
    
    # Count overlapping rectangles
    overlaps = find_overlaps()
    total_count = len(rectangles)
    
    # Add extra count for each overlap point that's shared between different rectangles
    for y, x in overlaps:
        overlap_count = 0
        for rect in rectangles:
            top, left, bottom, right = rect
            if (top <= y <= bottom and left <= x <= right):
                overlap_count += 1
        if overlap_count > 1:
            total_count += 1
    
    return total_count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

print(find_rectangles(grid))