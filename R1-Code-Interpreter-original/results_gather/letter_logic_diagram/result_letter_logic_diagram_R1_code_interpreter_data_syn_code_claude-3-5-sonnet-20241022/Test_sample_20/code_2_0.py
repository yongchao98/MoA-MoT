def get_diagonal_letter(grid):
    # Get the first letter found on the minor diagonal
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
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
    
    # If this is a diagonal cell, it MUST match the diagonal letter
    diagonal_letter = get_diagonal_letter(grid)
    if row + col == 6:  # Cell is on minor diagonal
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def find_empty(grid):
    # First priority: fill diagonal cells
    diagonal_letter = get_diagonal_letter(grid)
    if diagonal_letter:
        for i in range(7):
            j = 6-i
            if grid[i][j] == '':
                return i, j, True  # True indicates diagonal cell
    
    # Second priority: fill remaining cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j, False
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col, is_diagonal = empty
    letters = 'abcdefg'
    
    # If this is a diagonal cell or we have a diagonal letter, must use it
    diagonal_letter = get_diagonal_letter(grid)
    if is_diagonal and diagonal_letter:
        letters = diagonal_letter
    elif is_diagonal:
        # Try each letter for the first diagonal cell
        pass
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['b', '', '', 'g', '', 'c', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', '', 'a', 'c', '', 'b', 'f'],
    ['g', '', '', 'e', '', '', 'd'],
    ['', '', 'e', 'b', '', '', ''],
    ['', 'e', '', 'f', 'd', 'g', 'a'],
    ['e', 'b', 'f', 'd', 'g', '', 'c']
]

# First, identify the diagonal letter (e)
diagonal_letter = get_diagonal_letter(grid)
print(f"Diagonal letter must be: {diagonal_letter}")

if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")