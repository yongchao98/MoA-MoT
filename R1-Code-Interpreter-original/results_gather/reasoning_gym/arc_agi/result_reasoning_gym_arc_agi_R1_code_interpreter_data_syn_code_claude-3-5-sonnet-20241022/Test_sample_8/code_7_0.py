def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def is_enclosed_area(i, j):
        if grid[i][j] == 8:
            return False
            
        # Check if point is inside a U-shape or box of 8's
        
        # Check vertical U-shape (open to right)
        def check_vertical_u(x, y):
            if y > 0 and grid[x][y-1] == 8:  # Left wall
                if (x > 0 and grid[x-1][y-1] == 8 and  # Top-left corner
                    x < rows-1 and grid[x+1][y-1] == 8):  # Bottom-left corner
                    return True
            return False
            
        # Check horizontal U-shape (open to bottom)
        def check_horizontal_u(x, y):
            if x > 0 and grid[x-1][y] == 8:  # Top wall
                if (y > 0 and grid[x-1][y-1] == 8 and  # Top-left corner
                    y < cols-1 and grid[x-1][y+1] == 8):  # Top-right corner
                    return True
            return False
        
        # Check box formation
        def check_box(x, y):
            # Look for left and right walls
            left_wall = False
            right_wall = False
            for dy in [-1, 1]:
                ny = y + dy
                if 0 <= ny < cols and grid[x][ny] == 8:
                    if dy == -1:
                        left_wall = True
                    else:
                        right_wall = True
            
            # Look for top and bottom walls
            top_wall = False
            bottom_wall = False
            for dx in [-1, 1]:
                nx = x + dx
                if 0 <= nx < rows and grid[nx][y] == 8:
                    if dx == -1:
                        top_wall = True
                    else:
                        bottom_wall = True
            
            return (left_wall and right_wall) or (top_wall and bottom_wall)
        
        return check_vertical_u(i, j) or check_horizontal_u(i, j) or check_box(i, j)
    
    # Mark cells that should be filled
    for i in range(rows):
        for j in range(cols):
            if is_enclosed_area(i, j):
                output[i][j] = 0
    
    # Convert to string format
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,8,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,8,8,8,5,5,5,5],
    [5,8,8,8,8,5,5,5,5,5,5,5,5,5,5],
    [5,8,5,5,8,5,5,5,8,8,8,8,5,5,5],
    [5,8,5,5,8,5,5,5,8,5,5,8,5,5,5],
    [5,8,5,8,8,5,5,5,8,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,8,8,8,8,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
]

find_enclosed_areas(test_input)