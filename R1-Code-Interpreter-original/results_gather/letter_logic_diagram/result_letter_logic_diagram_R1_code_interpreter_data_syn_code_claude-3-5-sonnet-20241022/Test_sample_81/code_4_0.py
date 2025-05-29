def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid(grid, row, col, letter, initial_grid):
    # STRICT pre-filled check
    if initial_grid[row][col] != '':
        return letter == initial_grid[row][col]
    
    # Row check
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Column check
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Minor diagonal check (top-right to bottom-left)
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def solve(grid, initial_grid, pos=0):
    if pos >= 49:
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If this is a pre-filled position, must use that value
    if initial_grid[row][col] != '':
        if is_valid(grid, row, col, initial_grid[row][col], initial_grid):
            grid[row][col] = initial_grid[row][col]
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    # If on minor diagonal, must match existing diagonal value if any
    if row + col == 6:
        diagonal_value = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_value = grid[i][j]
                break
        if diagonal_value:
            if is_valid(grid, row, col, diagonal_value, initial_grid):
                grid[row][col] = diagonal_value
                if solve(grid, initial_grid, pos + 1):
                    return True
                grid[row][col] = ''
            return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, initial_grid):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid with pre-filled values
initial_grid = [
    ['', 'b', 'd', '', 'e', 'g', 'f'],
    ['', '', 'a', 'e', 'g', '', ''],
    ['d', 'a', '', '', '', 'c', ''],
    ['a', 'e', '', '', '', 'b', ''],
    ['', '', 'f', 'c', 'b', '', 'a'],
    ['', 'f', 'c', '', 'd', '', ''],
    ['f', '', 'b', '', '', '', 'g']
]

# Create empty working grid
grid = [['' for _ in range(7)] for _ in range(7)]

# First, verify all pre-filled positions
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            # Check if pre-filled positions violate row/column constraints
            for k in range(7):
                if k != j and initial_grid[i][k] == initial_grid[i][j]:
                    valid = False
                if k != i and initial_grid[k][j] == initial_grid[i][j]:
                    valid = False
            # Check if pre-filled positions on minor diagonal match
            if i + j == 6:
                for x in range(7):
                    y = 6 - x
                    if initial_grid[x][y] != '' and initial_grid[x][y] != initial_grid[i][j]:
                        valid = False

if not valid:
    print("Initial configuration is invalid!")
else:
    if solve(grid, initial_grid):
        print_solution(grid)
    else:
        print("No solution exists!")