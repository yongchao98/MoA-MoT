def check_initial_state(grid):
    # Check what letters are pre-filled on the diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    return diagonal_letters

def is_valid_for_diagonal(grid, letter):
    # Check if letter can be used on diagonal without violating constraints
    test_grid = [row[:] for row in grid]
    
    # Try placing the letter on all diagonal positions
    for i in range(7):
        j = 6-i
        # If position is filled and different, invalid
        if test_grid[i][j] != '' and test_grid[i][j] != letter:
            return False
        # If empty, check if letter appears elsewhere in row/column
        if test_grid[i][j] == '':
            if letter in test_grid[i] or letter in [test_grid[k][j] for k in range(7)]:
                return False
    return True

def fill_grid(grid, diagonal_letter):
    # First fill all diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    def get_next_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):  # Skip diagonal
                    return (i, j)
        return None
    
    def get_valid_letters(row, col):
        used_row = set(grid[row])
        used_col = set(grid[i][col] for i in range(7))
        return [c for c in 'abcdefg' if c not in used_row and c not in used_col]
    
    def solve():
        pos = get_next_empty()
        if not pos:
            return True
        
        row, col = pos
        for letter in get_valid_letters(row, col):
            if letter != diagonal_letter:  # Don't use diagonal letter elsewhere
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    return solve()

def solve_puzzle():
    initial = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    # Check pre-filled diagonal positions
    diagonal_letters = check_initial_state(initial)
    
    # If we have diagonal letters, they must be consistent
    if len(diagonal_letters) > 1:
        return None
    
    # Try each possible letter for diagonal
    for diagonal_letter in 'abcdefg':
        # If we have diagonal constraints, only try that letter
        if diagonal_letters and diagonal_letter not in diagonal_letters:
            continue
            
        if not is_valid_for_diagonal(initial, diagonal_letter):
            continue
            
        # Create a fresh grid and try to solve
        grid = [row[:] for row in initial]
        if fill_grid(grid, diagonal_letter):
            # Verify solution
            result = ""
            for row in grid:
                result += ','.join(row) + "\n"
            print(result.strip())
            return True
            
    return False

solve_puzzle()