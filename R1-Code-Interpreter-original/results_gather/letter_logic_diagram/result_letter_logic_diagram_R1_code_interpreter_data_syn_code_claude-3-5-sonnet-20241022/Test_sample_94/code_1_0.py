def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal - must all be the same letter
    if row + col == 6:  # If on minor diagonal
        # Find the required letter for minor diagonal
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
        if not diagonal_letter and row > 0:  # If first time filling diagonal
            # Check if this letter can fill the entire diagonal
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    return False
    return True

def find_empty(grid):
    # First fill minor diagonal positions
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, only try the letter that's already on the diagonal
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter:
            letters = diagonal_letter
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['a','f','g','c','e','d','b'],
    ['f','','c','','','b','a'],
    ['g','c','e','','','a','f'],
    ['c','','','','','f',''],
    ['e','d','','a','f','',''],
    ['d','b','a','','g','','e'],
    ['b','','','','','','d']
]

# Solve the puzzle
if solve_grid(initial_grid):
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")