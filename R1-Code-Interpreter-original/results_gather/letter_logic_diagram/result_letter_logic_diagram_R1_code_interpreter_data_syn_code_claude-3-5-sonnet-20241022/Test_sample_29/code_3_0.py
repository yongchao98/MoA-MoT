def get_possible_diagonal_letters(grid):
    # Check what letters are already on the minor diagonal
    diagonal_letters = set('abcdefg')
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            diagonal_letters = {grid[i][j]}
            break
    
    # Check what letters are still possible based on row/column constraints
    for i in range(7):
        j = 6 - i
        # Check row constraints
        for col in range(7):
            if col != j and grid[i][col] != '':
                if grid[i][col] in diagonal_letters:
                    diagonal_letters.discard(grid[i][col])
        # Check column constraints
        for row in range(7):
            if row != i and grid[row][j] != '':
                if grid[row][j] in diagonal_letters:
                    diagonal_letters.discard(grid[row][j])
    
    return diagonal_letters

def is_valid_move(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, must match diagonal_letter
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
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return i, j
    # Then other positions
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
    if row + col == 6:
        # On diagonal - must use diagonal_letter
        if is_valid_move(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    else:
        # Not on diagonal - try all letters except diagonal_letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve_grid(grid, diagonal_letter):
                    return True
                grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

# Find possible diagonal letters
possible_diagonals = get_possible_diagonal_letters(grid)
solved = False

for diagonal_letter in possible_diagonals:
    test_grid = [row[:] for row in grid]
    if solve_grid(test_grid, diagonal_letter):
        for row in test_grid:
            print(','.join(row))
        solved = True
        break

if not solved:
    print("No solution exists")