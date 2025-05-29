def quick_rectangle_count(grid):
    # Count rectangles by finding unique top-left corners
    count = 0
    rows = len(grid)
    cols = len(grid[0])
    
    # Count basic rectangles
    for i in range(rows):
        for j in range(cols):
            # Check if current position is a top-left corner
            if grid[i][j] in '#█':
                # Check if it's not part of a border (no border char above or left)
                is_corner = True
                if i > 0 and grid[i-1][j] in '#█':
                    is_corner = False
                if j > 0 and grid[i][j-1] in '#█':
                    is_corner = False
                
                if is_corner:
                    count += 1
    
    # Count overlapping sections (each █ represents one additional rectangle)
    overlap_count = sum(row.count('█') for row in grid)
    
    return count + overlap_count - 1  # Subtract 1 to avoid double counting

# Read input (limiting to relevant portion)
grid = []
start_row = None
end_row = None
start_col = None
end_col = None

# Read and find boundaries
row_num = 0
while True:
    try:
        line = input()
        if '#' in line or '█' in line:
            if start_row is None:
                start_row = row_num
            end_row = row_num
            
            # Find column boundaries
            for col, char in enumerate(line):
                if char in '#█':
                    if start_col is None:
                        start_col = col
                    start_col = min(start_col, col)
                    if end_col is None:
                        end_col = col
                    end_col = max(end_col, col)
        
        grid.append(line)
        row_num += 1
    except EOFError:
        break

# Trim grid to relevant portion
if start_row is not None and end_row is not None and start_col is not None and end_col is not None:
    grid = [row[start_col:end_col+1] for row in grid[start_row:end_row+1]]

print(f"<<<{quick_rectangle_count(grid)}>>>")