def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    # Compress grid into a binary representation for faster access
    binary_grid = []
    for row in grid:
        binary_row = []
        for char in row:
            binary_row.append(1 if char in ['#', '█'] else 0)
        binary_grid.append(binary_row)
    
    count = 0
    # Store last seen border position for each column
    last_border = [-1] * width
    
    # Single pass through the grid
    for y in range(height):
        start_x = -1
        for x in range(width):
            if binary_grid[y][x]:  # Found a border point
                if start_x == -1:
                    start_x = x
                last_border[x] = y
            elif start_x != -1:  # End of horizontal line
                # Check if we can form a rectangle with previous border points
                min_last = min(last_border[i] for i in range(start_x, x))
                if min_last != -1 and min_last < y:
                    # Verify if it forms a valid rectangle
                    valid = True
                    # Check bottom edge
                    for i in range(start_x, x):
                        if not binary_grid[min_last][i]:
                            valid = False
                            break
                    if valid:
                        count += 1
                start_x = -1
        
        # Handle case where line extends to end of row
        if start_x != -1:
            min_last = min(last_border[i] for i in range(start_x, width))
            if min_last != -1 and min_last < y:
                valid = True
                for i in range(start_x, width):
                    if not binary_grid[min_last][i]:
                        valid = False
                        break
                if valid:
                    count += 1
    
    return count

# Read grid efficiently
grid = []
try:
    while True:
        line = input()
        if line.strip():
            grid.append(line)
except EOFError:
    pass

# Count rectangles
print(f"<<<{find_rectangles(grid)}>>>")