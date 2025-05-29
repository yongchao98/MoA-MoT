def fill_diagonal_g(grid):
    # Fill all minor diagonal positions with 'g'
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = 'g'
        elif grid[i][j] != 'g':
            return False
    return True

def is_valid(grid, row, col, letter):
    # Skip validation for pre-filled cells
    if grid[row][col] != '':
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
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_grid(grid):
    # First, fill all diagonal positions with 'g'
    if not fill_diagonal_g(grid):
        return False
    
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for letter in 'abcdefg':
        if row + col == 6:  # Skip diagonal positions
            continue
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', 'g', '', ''],
    ['e', '', 'c', '', 'b', 'f', ''],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

# Create a copy of the initial grid to preserve pre-filled values
initial_grid = [[cell for cell in row] for row in grid]

# First verify that 'g' can be placed on all diagonal positions
test_grid = [[cell for cell in row] for row in initial_grid]
if fill_diagonal_g(test_grid):
    # Now solve with the diagonal constraint enforced
    if solve_grid(grid):
        result = ""
        for row in grid:
            result += ','.join(row) + "\n"
        print("<<<")
        print(result.strip())
        print(">>>")
    else:
        print("No solution exists")
else:
    print("Cannot place 'g' on diagonal due to conflicts")