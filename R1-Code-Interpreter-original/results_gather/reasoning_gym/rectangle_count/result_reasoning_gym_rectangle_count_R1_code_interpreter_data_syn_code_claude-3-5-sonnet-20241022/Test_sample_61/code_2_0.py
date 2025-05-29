def count_rectangles_fast(grid):
    # Quick count of rectangles and overlaps
    basic_count = 0
    overlap_count = 0
    
    # Count basic rectangles by identifying top-left corners
    # and overlap points in a single pass
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            # Count overlap points
            if grid[i][j] == '█':
                overlap_count += 1
            
            # Identify top-left corners of rectangles
            if grid[i][j] in '#█':
                # Check if it's a top-left corner
                is_corner = (i == 0 or grid[i-1][j] not in '#█') and (j == 0 or grid[i][j-1] not in '#█')
                if is_corner:
                    basic_count += 1

    # Each overlap point indicates an additional rectangle
    total = basic_count + overlap_count

    print(total)

# Read input efficiently
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines at the start
            continue
        if '#' in line or '█' in line:  # Only store lines with relevant characters
            grid.append(line)
except EOFError:
    pass

count_rectangles_fast(grid)