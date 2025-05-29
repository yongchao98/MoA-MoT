def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial_values(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    return True

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

def solve_puzzle():
    # Initialize empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # Initial grid with given values
    initial_grid = [
        ['', '', '', '', 'd', 'g', 'c'],
        ['', 'a', '', 'd', '', 'c', ''],
        ['a', '', 'd', 'g', '', '', ''],
        ['', '', 'g', 'c', '', '', 'a'],
        ['d', '', 'c', 'e', '', '', 'b'],
        ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
        ['c', 'e', 'f', 'a', '', 'd', 'g']
    ]
    
    # First, fill the minor diagonal with 'c'
    for i in range(7):
        grid[i][6-i] = 'c'
    
    def solve(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return check_initial_values(grid, initial_grid)
        
        # Skip if cell is on minor diagonal or pre-filled
        if col == 6-row or (initial_grid[row][col] != '' and grid[row][col] != ''):
            return solve(row, col + 1)
        
        # If position is pre-filled
        if initial_grid[row][col] != '':
            letter = initial_grid[row][col]
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
            grid[row][col] = ''
            return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if letter != 'c' and is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve():
        return grid
    return None

# Solve and print the puzzle
result = solve_puzzle()
if result:
    print_grid(result)
else:
    print("No solution found")