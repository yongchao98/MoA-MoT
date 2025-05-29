def is_valid_grid(grid):
    size = 7
    letters = set('abcdefg')
    
    # Check rows
    for row in grid:
        if row.count('') == 0 and set(row) != letters:
            return False
    
    # Check columns
    for j in range(size):
        col = [grid[i][j] for i in range(size)]
        if col.count('') == 0 and set(col) != letters:
            return False
    
    return True

def get_next_empty(grid, size=7):
    # First fill non-diagonal cells
    for i in range(size):
        for j in range(size):
            if grid[i][j] == '' and (i + j != size - 1):
                return (i, j)
    return None

def get_valid_letters(grid, row, col, size=7):
    used_row = set(x for x in grid[row] if x != '')
    used_col = set(grid[i][col] for i in range(size) if grid[i][col] != '')
    return set('abcdefg') - used_row - used_col

def solve_puzzle():
    size = 7
    # Initialize grid with given values
    grid = [
        ['a', 'f', 'b', 'd', '', 'c', 'e'],
        ['', '', 'd', '', 'c', '', 'e'],
        ['', '', 'e', '', '', '', 'e'],
        ['', '', 'c', '', '', 'f', 'e'],
        ['', '', '', 'a', '', 'b', 'e'],
        ['', '', '', 'f', 'b', '', 'e'],
        ['', 'a', 'f', '', '', 'e', 'e']
    ]
    
    def solve(grid):
        if not is_valid_grid(grid):
            return False
            
        pos = get_next_empty(grid)
        if not pos:
            return True
            
        row, col = pos
        valid_letters = get_valid_letters(grid, row, col)
        
        for letter in sorted(valid_letters):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
            
        return False
    
    # First, fill minor diagonal with 'e' (trying a different letter)
    for i in range(size):
        if grid[i][size-1-i] == '':
            grid[i][size-1-i] = 'e'
    
    if solve(grid):
        print("<<<")
        for row in grid:
            print(','.join(row))
        print(">>>")
        return True
    return False

# Try to solve
solve_puzzle()