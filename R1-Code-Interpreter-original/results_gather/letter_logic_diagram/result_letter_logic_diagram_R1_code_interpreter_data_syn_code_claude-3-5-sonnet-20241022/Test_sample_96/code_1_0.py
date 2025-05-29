def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check if it's on minor diagonal
    if row + col == 6:
        # First, find what letter should be on the diagonal
        diag_letter = None
        for i in range(7):
            for j in range(7):
                if i + j == 6 and grid[i][j] != '':
                    diag_letter = grid[i][j]
                    break
            if diag_letter:
                break
        
        # If we found a diagonal letter, enforce it
        if diag_letter and letter != diag_letter:
            return False
        # If we haven't found a diagonal letter yet, prefer 'e'
        if not diag_letter and letter != 'e':
            return False
    
    return True

def solve_grid(grid, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return True
    
    if grid[row][col] != '':
        return solve_grid(grid, row, col + 1)
    
    letters = 'abcdefg'
    # If on diagonal, try the diagonal letter first
    if row + col == 6:
        diag_letter = None
        for i in range(7):
            for j in range(7):
                if i + j == 6 and grid[i][j] != '':
                    diag_letter = grid[i][j]
                    break
            if diag_letter:
                break
        if diag_letter:
            letters = diag_letter + ''.join(l for l in 'abcdefg' if l != diag_letter)
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['g','','','','','',''],
    ['','','c','','','b','g'],
    ['d','','','','','g',''],
    ['c','','','','','f',''],
    ['e','','b','g','f','','c'],
    ['a','b','','f','','',''],
    ['','','f','','c','','a']
]

# Solve the puzzle
if solve_grid(initial_grid):
    for row in initial_grid:
        print(','.join(row))
else:
    print("No solution found")