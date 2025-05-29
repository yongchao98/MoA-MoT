def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check if letter already exists in row
    for j in range(7):
        if grid[row][j] == letter:
            return False
            
    # Check if letter already exists in column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    # Check minor diagonal constraint
    if row + col == 6:  # if on minor diagonal
        return letter == 'd'
        
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        if grid[i][6-i] == '':
            return (i, 6-i)
    
    # Then fill rest of grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        return True
        
    row, col = empty
    
    # If on minor diagonal, only try 'd'
    if row + col == 6:
        if is_valid(grid, row, col, 'd'):
            grid[row][col] = 'd'
            if solve(grid):
                return True
            grid[row][col] = ''
    else:
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# Verify pre-filled cells don't violate constraints
valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            temp = grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, temp):
                valid = False
                break
            grid[i][j] = temp
    if not valid:
        break

if valid and solve(grid):
    print_grid(grid)
else:
    print("No solution found")