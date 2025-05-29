def is_valid_move(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, letter must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
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

def solve_grid(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, only try the diagonal_letter
    if row + col == 6:
        letters = diagonal_letter
    
    for letter in letters:
        if is_valid_move(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

def try_solution(initial_grid, diagonal_letter):
    # Create a copy of the initial grid
    grid = [row[:] for row in initial_grid]
    
    # First verify if diagonal_letter is compatible with existing cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != diagonal_letter:
            return None
        grid[i][j] = diagonal_letter
    
    # Try to solve with this diagonal letter
    if solve_grid(grid, diagonal_letter):
        return grid
    return None

# Initial grid
initial_grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

# Try each possible letter for the diagonal
for diagonal_letter in 'abcdefg':
    solution = try_solution(initial_grid, diagonal_letter)
    if solution:
        for row in solution:
            print(','.join(row))
        break