def print_grid(grid):
    for row in grid:
        print(','.join(row))

ORIGINAL = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

def verify_solution(grid):
    # Verify all original values are preserved
    for i in range(7):
        for j in range(7):
            if ORIGINAL[i][j] != '' and grid[i][j] != ORIGINAL[i][j]:
                return False
    
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:
        return False
    
    return True

def solve():
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Copy original values
    for i in range(7):
        for j in range(7):
            if ORIGINAL[i][j] != '':
                grid[i][j] = ORIGINAL[i][j]
    
    # Find diagonal letter
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    def get_available_letters(row, col):
        used_row = set(grid[row])
        used_col = set(grid[i][col] for i in range(7))
        
        if row + col == 6:  # On minor diagonal
            if diag_letter:
                return {diag_letter}
            else:
                # Find a letter that works for the entire diagonal
                for letter in 'abcdefg':
                    valid = True
                    for i in range(7):
                        if grid[i][6-i] != '' and grid[i][6-i] != letter:
                            valid = False
                            break
                    if valid:
                        return {letter}
        
        return set('abcdefg') - used_row - used_col
    
    def solve_cell(pos):
        if pos >= 49:
            return verify_solution(grid)
        
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if ORIGINAL[row][col] != '':
            return solve_cell(pos + 1)
        
        available = get_available_letters(row, col)
        for letter in available:
            grid[row][col] = letter
            if solve_cell(pos + 1):
                return True
            grid[row][col] = ''
        
        return False
    
    if solve_cell(0):
        print_grid(grid)
        return True
    else:
        print("No solution exists")
        return False

solve()