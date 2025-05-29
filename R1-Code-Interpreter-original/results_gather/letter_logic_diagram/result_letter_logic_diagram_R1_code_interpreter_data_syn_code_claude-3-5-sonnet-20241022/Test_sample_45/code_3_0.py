def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        # All minor diagonal positions must contain the same letter
        for i in range(7):
            if grid[i][6-i] != '' and grid[i][6-i] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    # If this is a minor diagonal position, we need to find what letter should go there
    if row + col == 6:
        # Check if we already have a letter on the minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            if is_valid(grid, row, col, diag_letter):
                grid[row][col] = diag_letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
        else:
            # Try each letter for diagonal
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    if solve_grid(grid):
                        return True
                    grid[row][col] = ''
    else:
        # Regular position
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_grid(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid with given values
grid = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

if solve_grid(grid):
    print_grid(grid)
else:
    print("No solution found")