def get_valid_letters(grid, row, col, initial_grid):
    # If cell is pre-filled, only that letter is valid
    if initial_grid[row][col] != '':
        return {initial_grid[row][col]}
    
    # If cell is on minor diagonal, only 'd' is valid
    if row + col == 6:
        return {'d'}
    
    # Get used letters in row and column
    used = set()
    for i in range(7):
        if grid[row][i] != '':
            used.add(grid[row][i])
        if grid[i][col] != '':
            used.add(grid[i][col])
    
    return set('abcdefg') - used

def verify_grid(grid, initial_grid):
    # Check if solution matches initial grid
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

def solve():
    initial = [
        ['b','','','','g','e','d'],
        ['','','','','','d','b'],
        ['','','','e','','',''],
        ['','','','','b','','f'],
        ['','','d','','','f',''],
        ['e','d','','','f','',''],
        ['d','b','','f','','g','e']
    ]
    
    # Create working grid
    grid = [[''] * 7 for _ in range(7)]
    
    # First, copy all initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
    
    # Then fill minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
    
    def solve_cell(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return verify_grid(grid, initial)
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return solve_cell(pos + 1)
        
        valid_letters = get_valid_letters(grid, row, col, initial)
        
        for letter in valid_letters:
            grid[row][col] = letter
            if solve_cell(pos + 1):
                return True
            grid[row][col] = ''
        
        return False
    
    if solve_cell():
        return grid
    return None

# Solve and output
solution = solve()
if solution:
    print('<<<')
    for row in solution:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")