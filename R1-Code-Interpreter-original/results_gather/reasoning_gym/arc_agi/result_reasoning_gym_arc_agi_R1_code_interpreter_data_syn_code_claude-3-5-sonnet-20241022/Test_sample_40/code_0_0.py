def find_sequences(grid):
    height = len(grid)
    width = len(grid[0])
    new_grid = [row[:] for row in grid]
    
    # Function to check if a position is valid
    def is_valid(x, y):
        return 0 <= x < width and 0 <= y < height
    
    # Function to process a sequence
    def process_sequence(positions, direction):
        if len(positions) < 3:  # Need at least 3 cells for a valid sequence
            return
        
        # Find the 6 position
        six_pos = None
        for pos in positions:
            y, x = pos
            if grid[y][x] == 6:
                six_pos = (y, x)
                break
        
        if six_pos is None:
            return
            
        # Determine sequence direction from 6
        six_y, six_x = six_pos
        
        # Add markers (1 and 3)
        if direction == 'horizontal':
            # Check left side
            if is_valid(six_x - 1, six_y) and grid[six_y][six_x - 1] == 4:
                if is_valid(six_x - 2, six_y) and grid[six_y][six_x - 2] == 0:
                    new_grid[six_y][six_x - 2] = 3
            # Check right side
            if is_valid(six_x + 1, six_y) and grid[six_y][six_x + 1] == 4:
                if is_valid(six_x + 2, six_y) and grid[six_y][six_x + 2] == 0:
                    new_grid[six_y][six_x + 2] = 1
                    
        elif direction == 'vertical':
            # Check above
            if is_valid(six_x, six_y - 1) and grid[six_y - 1][six_x] == 4:
                if is_valid(six_x, six_y - 2) and grid[six_y - 2][six_x] == 0:
                    new_grid[six_y - 2][six_x] = 3
            # Check below
            if is_valid(six_x, six_y + 1) and grid[six_y + 1][six_x] == 4:
                if is_valid(six_x, six_y + 2) and grid[six_y + 2][six_x] == 0:
                    new_grid[six_y + 2][six_x] = 1
    
    # Find horizontal sequences
    for y in range(height):
        sequence = []
        for x in range(width):
            if grid[y][x] in [4, 6]:
                sequence.append((y, x))
            elif sequence:
                process_sequence(sequence, 'horizontal')
                sequence = []
        if sequence:
            process_sequence(sequence, 'horizontal')
    
    # Find vertical sequences
    for x in range(width):
        sequence = []
        for y in range(height):
            if grid[y][x] in [4, 6]:
                sequence.append((y, x))
            elif sequence:
                process_sequence(sequence, 'vertical')
                sequence = []
        if sequence:
            process_sequence(sequence, 'vertical')
    
    # Format the output
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