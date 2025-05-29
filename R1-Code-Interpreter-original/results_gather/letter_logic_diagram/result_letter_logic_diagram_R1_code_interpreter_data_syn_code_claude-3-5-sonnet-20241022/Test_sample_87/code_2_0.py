def get_diagonal_letter(grid):
    # Get the first non-empty letter on minor diagonal
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
    
    # Check minor diagonal requirement
    diagonal_letter = get_diagonal_letter(grid)
    if row + col == 6:  # Cell is on minor diagonal
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve(grid):
    # Find empty cell
    empty = None
    # First try to fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            empty = (i, 6-i)
            break
    # Then other cells
    if not empty:
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    empty = (i, j)
                    break
            if empty:
                break
    
    if not empty:
        return True
    
    row, col = empty
    diagonal_letter = get_diagonal_letter(grid)
    
    # If on minor diagonal and we have a diagonal letter, only try that
    if row + col == 6 and diagonal_letter:
        letters = [diagonal_letter]
    # If on minor diagonal and no diagonal letter yet, try each letter
    elif row + col == 6:
        letters = 'abcdefg'
    # For non-diagonal cells, try all letters
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

if solve(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")