def verify_initial(grid, row, col, letter, initial):
    return initial[row][col] == '' or initial[row][col] == letter

def verify_diagonal(grid, row, col, letter):
    if row + col == 6:  # on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    return True

def get_diagonal_letter(grid):
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def solve():
    # Initial configuration
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    
    # Working grid
    grid = [[''] * 7 for _ in range(7)]
    
    # Copy initial values and verify they don't conflict
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
    
    # Track used letters in each row and column
    used_in_row = [set() for _ in range(7)]
    used_in_col = [set() for _ in range(7)]
    
    # Initialize tracking sets with pre-filled values
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                used_in_row[i].add(grid[i][j])
                used_in_col[j].add(grid[i][j])
    
    def is_valid(row, col, letter):
        # Check initial configuration
        if not verify_initial(grid, row, col, letter, initial):
            return False
            
        # Check row and column
        if letter in used_in_row[row] or letter in used_in_col[col]:
            return False
            
        # Check diagonal
        if row + col == 6:
            diag_letter = get_diagonal_letter(grid)
            if diag_letter and diag_letter != letter:
                return False
                
        return True
    
    def solve_cell(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return solve_cell(pos + 1)
        
        # Get available letters
        available = set('abcdefg') - used_in_row[row] - used_in_col[col]
        
        # If on diagonal, must match diagonal letter if one exists
        if row + col == 6:
            diag_letter = get_diagonal_letter(grid)
            if diag_letter:
                if diag_letter not in available:
                    return False
                available = {diag_letter}
        
        for letter in available:
            if is_valid(row, col, letter):
                grid[row][col] = letter
                used_in_row[row].add(letter)
                used_in_col[col].add(letter)
                
                if solve_cell(pos + 1):
                    return True
                    
                grid[row][col] = ''
                used_in_row[row].remove(letter)
                used_in_col[col].remove(letter)
                
        return False
    
    if solve_cell():
        for row in grid:
            print(','.join(row))
        return True
    return False

solve()