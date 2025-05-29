def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6:  # If cell is on minor diagonal
        # Get the letter that should be on the diagonal
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    letters = 'abcdefg'
    
    # If on minor diagonal, we need to match existing diagonal letter
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
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['b', '', '', '', '', 'a', 'f'],
    ['', '', '', '', '', '', 'b'],
    ['', 'e', 'd', 'a', '', '', ''],
    ['e', '', '', 'f', 'b', 'g', ''],
    ['', 'a', '', 'b', 'g', '', 'e'],
    ['', '', 'b', 'g', '', '', 'd'],
    ['', '', '', '', '', 'd', '']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution exists")