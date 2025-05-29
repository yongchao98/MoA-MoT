def get_valid_letters(grid, row, col, initial_grid):
    # If cell is pre-filled, ONLY return that letter
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
    
    # Create working grid and copy initial values
    grid = [[''] * 7 for _ in range(7)]
    
    # First, strictly copy all initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
                # Verify if this initial placement is valid
                row_vals = [x for x in grid[i] if x != '']
                if len(row_vals) != len(set(row_vals)):
                    return None
                col_vals = [grid[k][j] for k in range(7) if grid[k][j] != '']
                if len(col_vals) != len(set(col_vals)):
                    return None
    
    # Fill minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] == '' and (i, 6-i) not in [(r, c) for r in range(7) for c in range(7) if initial[r][c] != '']:
            grid[i][6-i] = 'd'
    
    def solve_cell(pos=0):
        if pos == 49:
            # Verify rows and columns
            for i in range(7):
                if set(grid[i]) != set('abcdefg'):  # Check row
                    return False
                if set(grid[j][i] for j in range(7)) != set('abcdefg'):  # Check column
                    return False
            return True
            
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