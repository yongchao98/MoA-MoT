def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, diag_letter):
    # Check if position is on minor diagonal
    if row + col == 6:
        if letter != diag_letter:
            return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve():
    # Initialize grid with pre-filled values
    grid = [
        ['g','','','','','',''],
        ['','','c','','','b','g'],
        ['d','','','','','g',''],
        ['c','','','','','f',''],
        ['e','','b','g','f','','c'],
        ['a','b','','f','','',''],
        ['','','f','','c','','a']
    ]
    
    # Find diagonal letter from pre-filled cells
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    if not diag_letter:
        diag_letter = 'e'  # If no diagonal letter found, use 'e'

    def backtrack():
        empty = find_empty(grid)
        if not empty:
            return True
        
        row, col = empty
        
        # If on diagonal, only try diagonal letter
        if row + col == 6:
            if is_valid(grid, row, col, diag_letter, diag_letter):
                grid[row][col] = diag_letter
                if backtrack():
                    return True
                grid[row][col] = ''
        else:
            # Try all possible letters
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter, diag_letter):
                    grid[row][col] = letter
                    if backtrack():
                        return True
                    grid[row][col] = ''
        return False

    if backtrack():
        print_grid(grid)
        return True
    return False

# Run the solver
if not solve():
    print("No solution found")