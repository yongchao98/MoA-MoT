def is_valid_state(grid, row, col, initial_grid):
    # Check if we're violating initial grid values
    if initial_grid[row][col] != '' and grid[row][col] != initial_grid[row][col]:
        return False
    
    # Check row up to current column
    row_letters = [x for x in grid[row][:col+1] if x != '']
    if len(row_letters) != len(set(row_letters)):
        return False
        
    # Check column up to current row
    col_letters = [grid[i][col] for i in range(row+1) if grid[i][col] != '']
    if len(col_letters) != len(set(col_letters)):
        return False
    
    return True

def solve_grid():
    initial_grid = [
        ['', '', '', '', 'd', 'g', 'c'],
        ['', 'a', '', 'd', '', 'c', ''],
        ['a', '', 'd', 'g', '', '', ''],
        ['', '', 'g', 'c', '', '', 'a'],
        ['d', '', 'c', 'e', '', '', 'b'],
        ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
        ['c', 'e', 'f', 'a', '', 'd', 'g']
    ]
    
    # Start with empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # First, fill the minor diagonal with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'
    
    def fill_remaining(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip diagonal positions and pre-filled positions
        if row + col == 6 or (initial_grid[row][col] != '' and grid[row][col] != ''):
            return fill_remaining(pos + 1)
        
        # If position is pre-filled in initial grid
        if initial_grid[row][col] != '':
            grid[row][col] = initial_grid[row][col]
            if is_valid_state(grid, row, col, initial_grid):
                if fill_remaining(pos + 1):
                    return True
            grid[row][col] = ''
            return False
        
        # Try each letter except 'c' (reserved for diagonal)
        for letter in 'abdefg':
            grid[row][col] = letter
            if is_valid_state(grid, row, col, initial_grid):
                if fill_remaining(pos + 1):
                    return True
            grid[row][col] = ''
        
        return False
    
    if fill_remaining():
        return grid
    return None

# Solve and print
result = solve_grid()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution found")