def count_rectangles(grid):
    # Count overlapping points
    overlap_count = sum(row.count('█') for row in grid)
    
    # Basic rectangle count (looking at continuous sections)
    base_count = 0
    
    # Count distinct rectangles by looking at top-left corners
    height = len(grid)
    width = len(grid[0])
    
    for y in range(height):
        for x in range(width):
            # Check if this could be a top-left corner
            if grid[y][x] in '#█':
                # Verify it's a top-left corner
                is_corner = (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█')
                if is_corner:
                    base_count += 1

    # Each overlapping point (█) represents two rectangles sharing that point
    # So we need to add half of the overlap count to account for additional rectangles
    total_count = base_count + (overlap_count // 2)
    
    return total_count

# Read grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(line)
except EOFError:
    pass

# Get result
print(count_rectangles(grid))