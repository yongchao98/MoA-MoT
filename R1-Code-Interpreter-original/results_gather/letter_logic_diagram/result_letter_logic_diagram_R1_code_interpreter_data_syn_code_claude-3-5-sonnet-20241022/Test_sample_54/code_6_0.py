def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If on minor diagonal, must match other filled diagonal positions
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates diagonal position
    
    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)  # False indicates non-diagonal position
    
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    
    # For diagonal positions, try to find existing diagonal letter
    if is_diagonal:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter:
            if is_valid(grid, row, col, diagonal_letter):
                grid[row][col] = diagonal_letter
                if solve(grid):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', '', '', '', 'd', ''],
    ['', 'g', '', '', '', 'd', 'f'],
    ['', '', '', '', 'f', 'e', 'a'],
    ['b', '', '', '', 'f', 'e', ''],
    ['c', '', '', '', 'a', 'g', ''],
    ['d', 'f', '', '', 'g', '', ''],
    ['', 'e', 'a', '', 'b', '', '']
]

# Solve the puzzle
if solve(initial_grid):
    result = []
    for row in initial_grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")