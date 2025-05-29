def find_forced_diagonal_letter(grid):
    letters = set('abcdefg')
    # Check existing diagonal positions
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    
    # If no diagonal position is filled, analyze constraints
    possible_letters = set('abcdefg')
    for i in range(7):
        j = 6-i
        # Remove letters that appear in the same row
        row_letters = set(grid[i]) - {''}
        possible_letters -= row_letters
        # Remove letters that appear in the same column
        col_letters = set(grid[r][j] for r in range(7)) - {''}
        possible_letters -= col_letters
    
    if len(possible_letters) == 1:
        return list(possible_letters)[0]
    return None

def is_valid(grid):
    # Check rows
    for row in grid:
        if len(set(row) - {''}) != len([x for x in row if x != '']):
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if len(set(col) - {''}) != len([x for x in col if x != '']):
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    filled_diagonal = [x for x in diagonal if x != '']
    if len(set(filled_diagonal)) > 1:
        return False
    
    return True

def solve(grid, diagonal_letter):
    if not is_valid(grid):
        return False
    
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    i, j = empty
    # If this is a diagonal position, only try diagonal_letter
    if i + j == 6:
        if is_valid_move(grid, i, j, diagonal_letter):
            grid[i][j] = diagonal_letter
            if solve(grid, diagonal_letter):
                return True
            grid[i][j] = ''
        return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if letter != diagonal_letter or i + j != 6:  # Don't use diagonal_letter except on diagonal
            if is_valid_move(grid, i, j, letter):
                grid[i][j] = letter
                if solve(grid, diagonal_letter):
                    return True
                grid[i][j] = ''
    
    return False

def is_valid_move(grid, i, j, letter):
    # Check row
    if letter in grid[i]:
        return False
    
    # Check column
    if letter in [grid[r][j] for r in range(7)]:
        return False
    
    return True

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

# First, find what letter must be on the diagonal
diagonal_letter = find_forced_diagonal_letter(grid)

# Fill all diagonal positions with this letter
test_grid = [row[:] for row in grid]
for i in range(7):
    j = 6-i
    if test_grid[i][j] == '':
        test_grid[i][j] = diagonal_letter
    elif test_grid[i][j] != diagonal_letter:
        print("No solution exists")
        exit()

if solve(test_grid, diagonal_letter):
    for row in test_grid:
        print(','.join(row))
else:
    print("No solution exists")