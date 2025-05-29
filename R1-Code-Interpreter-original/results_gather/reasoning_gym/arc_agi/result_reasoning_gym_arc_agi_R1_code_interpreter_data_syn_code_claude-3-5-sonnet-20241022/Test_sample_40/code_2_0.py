def find_sequences(grid):
    height = len(grid)
    width = len(grid[0])
    new_grid = [row[:] for row in grid]
    
    def is_valid(x, y):
        return 0 <= x < width and 0 <= y < height
    
    def process_sequence(positions, direction):
        if len(positions) < 2:
            return
            
        # Find if there's a 6 and its position
        six_pos = None
        sequence_positions = []
        for pos in positions:
            y, x = pos
            if grid[y][x] == 6:
                six_pos = (y, x)
            if grid[y][x] == 4:
                sequence_positions.append((y, x))
                
        if not six_pos:
            return
            
        six_y, six_x = six_pos
        
        if direction == 'horizontal':
            # If 6 is on the left side
            if any(x > six_x for _, x in sequence_positions):
                # Add 1's to the right of 4's
                for y, x in sequence_positions:
                    if x > six_x and is_valid(x+1, y) and grid[y][x+1] == 0:
                        new_grid[y][x+1] = 1
            
            # If 6 is on the right side
            if any(x < six_x for _, x in sequence_positions):
                # Add 3 to the left of the leftmost 4
                left_x = min(x for _, x in sequence_positions)
                if is_valid(left_x-1, six_y) and grid[six_y][left_x-1] == 0:
                    new_grid[six_y][left_x-1] = 3
                
        elif direction == 'vertical':
            # If 6 is on the top
            if any(y > six_y for y, _ in sequence_positions):
                # Add 1's below 4's
                for y, x in sequence_positions:
                    if y > six_y and is_valid(x, y+1) and grid[y+1][x] == 0:
                        new_grid[y+1][x] = 1
            
            # If 6 is on the bottom
            if any(y < six_y for y, _ in sequence_positions):
                # Add 3 above the topmost 4
                top_y = min(y for y, _ in sequence_positions)
                if is_valid(six_x, top_y-1) and grid[top_y-1][six_x] == 0:
                    new_grid[top_y-1][six_x] = 3
    
    # Process horizontal sequences
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
    
    # Process vertical sequences
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