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
        for pos in positions:
            y, x = pos
            if grid[y][x] == 6:
                six_pos = (y, x)
                break
                
        if not six_pos:
            return
            
        six_y, six_x = six_pos
        
        if direction == 'horizontal':
            # Process sequence to the right of 6
            if six_x < width - 1:
                x = six_x + 1
                while x < width and grid[six_y][x] == 4:
                    if x + 1 < width and grid[six_y][x + 1] == 0:
                        new_grid[six_y][x + 1] = 1
                    x += 1
            
            # Process sequence to the left of 6
            if six_x > 0:
                x = six_x - 1
                while x >= 0 and grid[six_y][x] == 4:
                    if x - 1 >= 0 and grid[six_y][x - 1] == 0:
                        new_grid[six_y][x - 1] = 3
                    x -= 1
                    
        elif direction == 'vertical':
            # Process sequence below 6
            if six_y < height - 1:
                y = six_y + 1
                while y < height and grid[y][six_x] == 4:
                    if y + 1 < height and grid[y + 1][six_x] == 0:
                        new_grid[y + 1][six_x] = 1
                    y += 1
            
            # Process sequence above 6
            if six_y > 0:
                y = six_y - 1
                while y >= 0 and grid[y][six_x] == 4:
                    if y - 1 >= 0 and grid[y - 1][six_x] == 0:
                        new_grid[y - 1][six_x] = 3
                    y -= 1
    
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