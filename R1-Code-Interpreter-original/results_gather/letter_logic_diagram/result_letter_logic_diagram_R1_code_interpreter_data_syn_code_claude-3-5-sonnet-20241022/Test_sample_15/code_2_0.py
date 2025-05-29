def get_possible_diagonal_letters(grid):
    possible = set('abcdefg')
    # Check what letters are already used in rows and columns of diagonal positions
    for i in range(7):
        # Position on diagonal is (i, 6-i)
        # Check row i
        for j in range(7):
            if grid[i][j] != '' and j != 6-i:
                if grid[i][j] in possible:
                    possible.remove(grid[i][j])
        # Check column 6-i
        for j in range(7):
            if grid[j][6-i] != '' and j != i:
                if grid[j][6-i] in possible:
                    possible.remove(grid[j][6-i])
    return possible

def is_valid(grid, row, col, letter):
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

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    return False

# Initial grid
grid = [
    ['c','','','','a','','b'],
    ['g','e','f','a','','',''],
    ['','','','d','','','g'],
    ['f','a','','','c','g',''],
    ['','d','','c','','e','f'],
    ['d','','','','','f',''],
    ['','','','','','','']
]

# First, find possible letters for diagonal
possible_diagonal = get_possible_diagonal_letters(grid)
print("Possible diagonal letters:", possible_diagonal)

for diag_letter in possible_diagonal:
    # Create a copy of the grid
    test_grid = [row[:] for row in grid]
    
    # Fill diagonal with the letter
    for i in range(7):
        if test_grid[i][6-i] == '':
            test_grid[i][6-i] = diag_letter
    
    # Try to solve the rest
    if solve(test_grid):
        print("<<<")
        for row in test_grid:
            print(','.join(row))
        print(">>>")
        break