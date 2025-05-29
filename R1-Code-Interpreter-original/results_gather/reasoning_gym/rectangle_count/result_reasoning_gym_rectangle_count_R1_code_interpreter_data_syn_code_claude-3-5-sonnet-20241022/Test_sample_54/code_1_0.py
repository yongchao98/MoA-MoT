def count_rectangles(grid):
    # Convert the grid into a list of strings
    grid = grid.split('\n')
    
    height = len(grid)
    if height == 0:
        return 0
    width = len(grid[0])
    
    count = 0
    
    # For each possible top-left corner
    for y1 in range(height-1):
        for x1 in range(width-1):
            if grid[y1][x1] != '#':
                continue
                
            # For each possible bottom-right corner
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if grid[y2][x2] != '#':
                        continue
                        
                    # Check if this forms a valid rectangle
                    is_rectangle = True
                    
                    # Check top and bottom edges
                    for x in range(x1, x2 + 1):
                        if grid[y1][x] != '#' or grid[y2][x] != '#':
                            is_rectangle = False
                            break
                            
                    if not is_rectangle:
                        continue
                        
                    # Check left and right edges
                    for y in range(y1, y2 + 1):
                        if grid[y][x1] != '#' or grid[y][x2] != '#':
                            is_rectangle = False
                            break
                    
                    if is_rectangle:
                        count += 1

    return count // 4  # Divide by 4 as each rectangle is counted from each corner

# Create test grid (80x80 with the rectangle)
grid = [" " * 80] * 40  # First 40 empty lines
rectangle_start = [" " * 55 + "#" * 9] * 1  # Top of rectangle
rectangle_middle = [" " * 55 + "#" + " " * 7 + "#"] * 1  # Middle of rectangle
rectangle_end = [" " * 55 + "#" * 9] * 1  # Bottom of rectangle
grid = grid + rectangle_start + rectangle_middle + rectangle_end + [" " * 80] * 37  # Rest of empty lines

# Join the grid lines
grid = "\n".join(grid)

print(count_rectangles(grid))