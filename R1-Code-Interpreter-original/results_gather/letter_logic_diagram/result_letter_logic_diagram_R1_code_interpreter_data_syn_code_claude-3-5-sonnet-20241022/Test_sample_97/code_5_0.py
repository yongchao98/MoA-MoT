def print_grid(grid):
    for row in grid:
        print(','.join(row))

# Initial configuration as a global constant
INITIAL = [
    ['', 'a', 'e', 'g', 'f', 'b', ''],
    ['a', '', 'g', 'f', '', '', ''],
    ['', '', '', '', 'c', 'd', 'a'],
    ['', 'f', 'b', 'c', '', 'a', ''],
    ['', 'b', 'c', '', 'a', '', 'g'],
    ['b', '', '', '', '', 'g', ''],
    ['c', '', '', 'e', 'g', '', '']
]

def verify_solution(grid):
    # Verify initial configuration
    for i in range(7):
        for j in range(7):
            if INITIAL[i][j] != '' and grid[i][j] != INITIAL[i][j]:
                return False
    
    # Verify rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[r][i] for r in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Verify diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def solve():
    grid = [[''] * 7 for _ in range(7)]
    
    # First, fill in all initial values
    for i in range(7):
        for j in range(7):
            if INITIAL[i][j] != '':
                grid[i][j] = INITIAL[i][j]
    
    def get_used_letters(row, col):
        # Get letters used in row
        row_letters = set(grid[row])
        # Get letters used in column
        col_letters = set(grid[i][col] for i in range(7))
        # Get diagonal letter if applicable
        diag_letter = None
        if row + col == 6:
            for i in range(7):
                if grid[i][6-i] != '':
                    diag_letter = grid[i][6-i]
                    break
        return row_letters, col_letters, diag_letter
    
    def is_valid_move(row, col, letter):
        if INITIAL[row][col] != '' and INITIAL[row][col] != letter:
            return False
        
        row_letters, col_letters, diag_letter = get_used_letters(row, col)
        if letter in row_letters or letter in col_letters:
            return False
        
        if row + col == 6 and diag_letter and letter != diag_letter:
            return False
        
        return True
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def solve_recursive():
        empty = find_empty()
        if not empty:
            return verify_solution(grid)
        
        row, col = empty
        row_letters, col_letters, diag_letter = get_used_letters(row, col)
        
        # If on diagonal and we have a diagonal letter, must use it
        if row + col == 6 and diag_letter:
            if is_valid_move(row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve_recursive():
                    return True
                grid[row][col] = ''
            return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_move(row, col, letter):
                grid[row][col] = letter
                if solve_recursive():
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve_recursive():
        print_grid(grid)
        return True
    return False

solve()