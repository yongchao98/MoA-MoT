def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, original):
    # Must match original if cell was pre-filled
    if original[row][col] != '' and letter != original[row][col]:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter and j != col:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter and i != row:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        # Find first non-empty cell in minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_puzzle():
    # Original puzzle
    original = [
        ['a','e','','f','g','c',''],
        ['','b','','','','d','a'],
        ['b','f','g','c','d','a','e'],
        ['f','','c','','a','','b'],
        ['g','c','d','a','e','b','f'],
        ['c','','','e','b','','g'],
        ['d','','','','f','g','c']
    ]
    
    # Working grid
    grid = [row[:] for row in original]
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if original[row][col] != '':
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, original):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if backtrack():
        print_grid(grid)
        return True
    return False

solve_puzzle()