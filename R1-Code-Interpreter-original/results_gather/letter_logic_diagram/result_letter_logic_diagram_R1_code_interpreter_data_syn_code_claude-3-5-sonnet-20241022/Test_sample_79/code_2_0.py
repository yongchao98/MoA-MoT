def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_diagonal_letter(grid):
    # Get the letter that appears most on the minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    counts = {}
    for letter in diagonal:
        if letter != '':
            counts[letter] = counts.get(letter, 0) + 1
    if counts:
        return max(counts.items(), key=lambda x: x[1])[0]
    return None

def get_empty_cell(grid):
    # First fill minor diagonal cells if any are empty
    for i in range(7):
        if grid[i][6-i] == '':
            return i, 6-i
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
    return None

def is_valid(grid, row, col, letter, diag_letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        return letter == diag_letter
    
    return True

def solve(grid, diag_letter):
    empty = get_empty_cell(grid)
    if not empty:
        return True
    
    row, col = empty
    
    # If on minor diagonal, only try the diagonal letter
    if row + col == 6:
        letters = [diag_letter]
    else:
        letters = 'abcdefg'
    
    for letter in letters:
        if is_valid(grid, row, col, letter, diag_letter):
            grid[row][col] = letter
            if solve(grid, diag_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['d','e','c','a','f','','g'],
    ['','c','','','b','g',''],
    ['c','a','f','','','d','e'],
    ['a','','','g','d','e','c'],
    ['f','','g','d','e','c','a'],
    ['b','g','','','c','','f'],
    ['g','d','e','c','','f','']
]

# First determine which letter should be on the minor diagonal
diag_letter = get_diagonal_letter(grid)
if diag_letter is None:
    diag_letter = 'b'  # We can see 'b' appears on the diagonal in the initial grid

if solve(grid, diag_letter):
    print_grid(grid)
else:
    print("No solution exists")