from itertools import permutations

# Initial grid setup
grid = [
    ['g', 'c', 'b', 'f', 'e', 'd', ''],
    ['', '', '', 'e', 'd', 'a', ''],
    ['', '', 'e', 'd', 'a', 'g', ''],
    ['f', 'e', '', '', 'g', 'c', 'b'],
    ['e', 'd', 'a', '', '', '', 'f'],
    ['d', '', 'g', 'c', '', '', 'e'],
    ['', 'g', '', 'b', 'f', '', '']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to solve the grid
def solve(grid):
    # Try each letter for the minor diagonal
    for diag_letter in 'abcdefg':
        # Create a copy of the grid to try this diagonal letter
        temp_grid = [row[:] for row in grid]
        # Place the diagonal letter
        for i in range(7):
            temp_grid[i][6-i] = diag_letter
        
        # Try to fill the rest of the grid
        if fill_grid(temp_grid, 0, 0):
            return temp_grid
    return None

# Recursive function to fill the grid
def fill_grid(grid, row, col):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row][col] != '':
        return fill_grid(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if can_place(grid, row, col, letter):
            grid[row][col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row][col] = ''
    return False

# Solve the grid
solution = solve(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")