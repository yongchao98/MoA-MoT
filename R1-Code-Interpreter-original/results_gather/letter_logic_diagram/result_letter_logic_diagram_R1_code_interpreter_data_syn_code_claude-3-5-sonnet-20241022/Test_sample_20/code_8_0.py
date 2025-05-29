def create_empty_grid():
    return [['' for _ in range(7)] for _ in range(7)]

def copy_initial_state(grid):
    initial = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]

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

def solve():
    # Start with empty grid
    grid = create_empty_grid()
    
    # First, fill all diagonal positions with 'e'
    for i in range(7):
        grid[i][6-i] = 'e'
    
    # Now try to place the initial values that don't conflict
    initial = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    
    # Copy non-diagonal initial values if they don't conflict
    for i in range(7):
        for j in range(7):
            if i + j != 6 and initial[i][j] != '':  # Skip diagonal positions
                if is_valid(grid, i, j, initial[i][j]):
                    grid[i][j] = initial[i][j]
                else:
                    return None  # Initial state is invalid
    
    def fill_remaining(row, col):
        if col >= 7:
            row += 1
            col = 0
        if row >= 7:
            return True
            
        # Skip filled positions and diagonal positions
        if grid[row][col] != '':
            return fill_remaining(row, col + 1)
            
        for letter in 'abcdefg':
            if row + col != 6 or letter == 'e':  # Only allow 'e' on diagonal
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    if fill_remaining(row, col + 1):
                        return True
                    grid[row][col] = ''
        return False
    
    if fill_remaining(0, 0):
        return grid
    return None

# Solve and print result
solution = solve()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")