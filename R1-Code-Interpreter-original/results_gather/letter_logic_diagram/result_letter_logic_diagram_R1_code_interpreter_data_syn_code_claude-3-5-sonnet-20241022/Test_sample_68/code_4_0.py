def get_diagonal_letter(grid):
    # Get the letter that must be on the diagonal (g from top-right corner)
    return grid[0][6]  # We know this is 'g' from the input

def is_valid(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, only diagonal_letter is allowed
    if row + col == 6:
        return letter == diagonal_letter
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def print_grid(grid):
    for row in grid:
        print(','.join(row))

def solve(grid, diagonal_letter):
    # Find empty cell
    empty_found = False
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                row, col = i, j
                empty_found = True
                break
        if empty_found:
            break
    
    if not empty_found:
        return True
    
    # Try filling the empty cell
    possible_letters = [diagonal_letter] if row + col == 6 else 'abcdefg'
    
    for letter in possible_letters:
        if is_valid(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', 'a', 'c', '', '', 'g'],
    ['', 'a', 'c', 'b', '', 'g', ''],
    ['a', '', '', '', 'g', '', ''],
    ['c', 'b', 'f', '', 'e', 'd', 'a'],
    ['b', 'f', '', '', '', 'a', ''],
    ['', '', 'e', '', '', 'c', ''],
    ['', '', 'd', 'a', '', 'b', 'f']
]

# First, fill all diagonal positions with 'g'
diagonal_letter = get_diagonal_letter(grid)
for i in range(7):
    j = 6 - i
    if grid[i][j] == '':
        grid[i][j] = diagonal_letter

# Now solve the rest of the grid
if solve(grid, diagonal_letter):
    print_grid(grid)
else:
    print("No solution exists")