def get_diagonal_letter(grid):
    # First, find what letters are already on the diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    
    # For each position on diagonal, eliminate impossible letters
    possible = set('abcdefg')
    for i in range(7):
        # Check row constraints
        for j in range(7):
            if j != 6-i and grid[i][j] != '':
                if grid[i][j] in possible:
                    possible.remove(grid[i][j])
        
        # Check column constraints
        j = 6-i
        for k in range(7):
            if k != i and grid[k][j] != '':
                if grid[k][j] in possible:
                    possible.remove(grid[k][j])
    
    # If there are pre-filled diagonal letters, they must match
    if diagonal_letters:
        possible &= diagonal_letters
    
    return list(possible)[0] if len(possible) > 0 else None

def is_valid(grid, row, col, letter, diagonal_letter):
    # If on diagonal, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve(grid, diagonal_letter):
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
    
    row, col = empty
    # If on diagonal, only try diagonal_letter
    if row + col == 6:
        if is_valid(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    else:
        # Try all letters except diagonal_letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve(grid, diagonal_letter):
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

# First, determine what letter must be on the diagonal
diagonal_letter = get_diagonal_letter(grid)
if diagonal_letter is None:
    print("No solution exists - cannot determine diagonal letter")
    exit()

print(f"Using diagonal letter: {diagonal_letter}")

# Fill all diagonal positions with this letter
test_grid = [row[:] for row in grid]
for i in range(7):
    if test_grid[i][6-i] == '':
        test_grid[i][6-i] = diagonal_letter

if solve(test_grid, diagonal_letter):
    print("<<<")
    for row in test_grid:
        print(','.join(row))
    print(">>>")
else:
    print("No solution exists")