def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_initial_grid():
    return [
        ['f', 'd', 'a', '.', '.', 'g', '.'],
        ['d', '.', '.', 'b', '.', '.', 'f'],
        ['.', 'c', 'b', '.', 'e', '.', '.'],
        ['.', '.', 'g', 'e', '.', '.', '.'],
        ['.', 'g', 'e', 'f', '.', 'a', 'c'],
        ['.', 'e', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'a', '.', '.', '.']
    ]

def get_available_values(grid, r, c, initial):
    used = set()
    
    # Check row
    for j in range(7):
        if grid[r][j] != '.':
            used.add(grid[r][j])
    
    # Check column
    for i in range(7):
        if grid[i][c] != '.':
            used.add(grid[i][c])
    
    # If on minor diagonal, must be 'a'
    if r + c == 6:
        return {'a'} - used
    
    return set('abcdefg') - used

def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '.' and grid[i][j] != initial[i][j]:
                return False
    return True

def verify_solution(grid):
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # Check row
            return False
        col = [grid[j][i] for j in range(7)]
        if len(set(col)) != 7:  # Check column
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == 'a' for x in diag):
        return False
    
    return True

def solve(grid, initial, r=0, c=0):
    if c == 7:
        r += 1
        c = 0
    if r == 7:
        return verify_solution(grid) and verify_initial(grid, initial)

    # Skip if cell is already filled
    if initial[r][c] != '.':
        grid[r][c] = initial[r][c]
        return solve(grid, initial, r, c + 1)

    available = get_available_values(grid, r, c, initial)
    for val in available:
        grid[r][c] = val
        if solve(grid, initial, r, c + 1):
            return True
        grid[r][c] = '.'
    
    return False

# Initialize grid
initial = get_initial_grid()
grid = [['.' for _ in range(7)] for _ in range(7)]

# Copy initial values
for i in range(7):
    for j in range(7):
        if initial[i][j] != '.':
            grid[i][j] = initial[i][j]

# First ensure all diagonal positions will be 'a'
for i in range(7):
    j = 6 - i
    if grid[i][j] != '.' and grid[i][j] != 'a':
        print("No solution exists - diagonal constraint cannot be satisfied")
        exit()

if solve(grid, initial):
    print_grid(grid)
else:
    print("No solution exists")