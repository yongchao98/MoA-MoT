def find_sequences(grid):
    height = len(grid)
    width = len(grid[0])
    new_grid = [row[:] for row in grid]
    
    def is_valid(x, y):
        return 0 <= x < width and 0 <= y < height
    
    def process_horizontal_sequence(y, start_x, end_x):
        # Find 6 position
        six_x = None
        for x in range(start_x, end_x + 1):
            if grid[y][x] == 6:
                six_x = x
                break
        
        if six_x is None:
            return
            
        # If 6 is on the right
        if six_x == end_x:
            # Add 3 at the left end
            if is_valid(start_x - 1, y) and grid[y][start_x - 1] == 0:
                new_grid[y][start_x - 1] = 3
            # Add 1's between 4's
            for x in range(start_x, end_x):
                if grid[y][x] == 4 and is_valid(x + 1, y) and grid[y][x + 1] == 0:
                    new_grid[y][x + 1] = 1
        
        # If 6 is on the left
        if six_x == start_x:
            # Add 3 at the right end
            if is_valid(end_x + 1, y) and grid[y][end_x + 1] == 0:
                new_grid[y][end_x + 1] = 3
            # Add 1's between 4's
            for x in range(start_x, end_x):
                if grid[y][x] == 4 and is_valid(x + 1, y) and grid[y][x + 1] == 0:
                    new_grid[y][x + 1] = 1
    
    def process_vertical_sequence(x, start_y, end_y):
        # Find 6 position
        six_y = None
        for y in range(start_y, end_y + 1):
            if grid[y][x] == 6:
                six_y = y
                break
        
        if six_y is None:
            return
            
        # If 6 is at the bottom
        if six_y == end_y:
            # Add 3 at the top end
            if is_valid(x, start_y - 1) and grid[start_y - 1][x] == 0:
                new_grid[start_y - 1][x] = 3
            # Add 1's between 4's
            for y in range(start_y, end_y):
                if grid[y][x] == 4 and is_valid(x, y + 1) and grid[y + 1][x] == 0:
                    new_grid[y + 1][x] = 1
        
        # If 6 is at the top
        if six_y == start_y:
            # Add 3 at the bottom end
            if is_valid(x, end_y + 1) and grid[end_y + 1][x] == 0:
                new_grid[end_y + 1][x] = 3
            # Add 1's between 4's
            for y in range(start_y, end_y):
                if grid[y][x] == 4 and is_valid(x, y + 1) and grid[y + 1][x] == 0:
                    new_grid[y + 1][x] = 1
    
    # Find horizontal sequences
    for y in range(height):
        start_x = None
        for x in range(width):
            if grid[y][x] in [4, 6]:
                if start_x is None:
                    start_x = x
            elif start_x is not None:
                process_horizontal_sequence(y, start_x, x - 1)
                start_x = None
        if start_x is not None:
            process_horizontal_sequence(y, start_x, width - 1)
    
    # Find vertical sequences
    for x in range(width):
        start_y = None
        for y in range(height):
            if grid[y][x] in [4, 6]:
                if start_y is None:
                    start_y = y
            elif start_y is not None:
                process_vertical_sequence(x, start_y, y - 1)
                start_y = None
        if start_y is not None:
            process_vertical_sequence(x, start_y, height - 1)
    
    result = '\n'.join(' '.join(str(cell) for cell in row) for row in new_grid)
    print(result)

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,1,3,0,0],
    [0,0,4,4,0,0,0,0,0,0,6,4,4,0,0],
    [0,0,4,0,4,0,0,0,0,1,4,0,4,0,0],
    [0,0,6,4,4,0,0,0,0,3,4,4,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,4,4,6,0,0,0,0,0,4,4,0],
    [0,0,0,0,4,0,4,0,0,0,0,4,0,4,0],
    [0,0,0,0,0,4,4,0,0,0,0,4,4,6,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

find_sequences(test_input)