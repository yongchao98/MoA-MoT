def verify_solution(grid):
    # Check rows and columns
    for i in range(7):
        row = set(grid[i])
        col = set(grid[j][i] for j in range(7))
        if len(row) != 7 or len(col) != 7:
            return False
    
    # Check diagonal
    diagonal = set(grid[i][6-i] for i in range(7))
    return len(diagonal) == 1

def create_initial_state():
    # Create grid with initial values
    grid = [
        ['f', 'd', '', '', '', 'c', 'g'],
        ['d', 'e', '', 'a', 'c', '', 'f'],
        ['e', 'b', '', '', 'g', 'f', 'g'],  # Force g in diagonal
        ['', '', 'c', 'g', 'f', 'g', 'e'],  # Force g in diagonal
        ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
        ['', 'g', '', 'd', 'e', 'g', 'a'],  # Force g in diagonal
        ['g', 'f', '', 'e', 'b', 'a', 'c']
    ]
    return grid

def is_valid_partial(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check diagonal requirement
    if row + col == 6 and letter != 'g':
        return False
    
    return True

def solve_puzzle(grid, pos=0):
    if pos == 49:
        return verify_solution(grid)
    
    row = pos // 7
    col = pos % 7
    
    # Skip pre-filled cells
    if grid[row][col] != '':
        return solve_puzzle(grid, pos + 1)
    
    # If on diagonal, only try 'g'
    if row + col == 6:
        if is_valid_partial(grid, row, col, 'g'):
            grid[row][col] = 'g'
            if solve_puzzle(grid, pos + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try each possible letter
    letters = 'abcdefg'
    for letter in letters:
        if is_valid_partial(grid, row, col, letter):
            grid[row][col] = letter
            if solve_puzzle(grid, pos + 1):
                return True
            grid[row][col] = ''
    
    return False

# Initialize and solve
grid = create_initial_state()

# First, ensure all diagonal positions are 'g'
for i in range(7):
    if grid[i][6-i] == '':
        grid[i][6-i] = 'g'
    elif grid[i][6-i] != 'g':
        print("Invalid initial state")
        exit()

if solve_puzzle(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")