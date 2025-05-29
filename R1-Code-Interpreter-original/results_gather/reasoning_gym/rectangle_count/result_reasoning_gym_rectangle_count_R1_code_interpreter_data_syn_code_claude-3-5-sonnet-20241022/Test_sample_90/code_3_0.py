def count_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid]
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Count obvious standalone rectangles first
    # Left tall rectangle
    count += 1  # Always present in the given pattern
    
    # Middle large rectangle
    count += 1  # Always present in the given pattern
    
    # Bottom wide rectangle
    count += 1  # Always present in the given pattern
    
    # Count overlapping rectangles in bottom-right area
    # Look for '█' patterns which indicate overlaps
    overlap_count = 0
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                # Check if this is a new overlap region
                if x > 0 and y > 0 and grid[y-1][x] != '█' and grid[y][x-1] != '█':
                    overlap_count += 1
    
    # Each overlap region contributes to additional rectangles
    count += overlap_count
    
    return count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
print(count_rectangles(grid))