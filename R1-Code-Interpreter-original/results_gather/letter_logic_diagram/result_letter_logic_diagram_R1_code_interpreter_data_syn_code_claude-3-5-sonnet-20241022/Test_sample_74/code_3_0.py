def get_diagonal_letter(initial_grid):
    # Find the first letter on the minor diagonal in the initial grid
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return 'c'  # We know it should be 'c' from the puzzle

def is_valid_placement(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col] != '']:
        return False
    
    # Check minor diagonal
    if row + col == 6 and letter != 'c':
        return False
        
    return True

def solve(grid, initial_grid, pos=0):
    if pos == 49:  # 7x7 = 49 cells
        return True
        
    row = pos // 7
    col = pos % 7
    
    # If cell is pre-filled, move to next cell
    if initial_grid[row][col] != '':
        grid[row][col] = initial_grid[row][col]
        return solve(grid, initial_grid, pos + 1)
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, initial_grid, pos + 1):
                return True
            grid[row][col] = ''
            
    return False

# Initialize grid with given values
initial_grid = [
    ['', '', '', '', 'd', 'g', 'c'],
    ['', 'a', '', 'd', '', 'c', ''],
    ['a', '', 'd', 'g', '', '', ''],
    ['', '', 'g', 'c', '', '', 'a'],
    ['d', '', 'c', 'e', '', '', 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', '', 'd', 'g']
]

# Create empty grid
grid = [[''] * 7 for _ in range(7)]

# Solve the puzzle
if solve(grid, initial_grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")