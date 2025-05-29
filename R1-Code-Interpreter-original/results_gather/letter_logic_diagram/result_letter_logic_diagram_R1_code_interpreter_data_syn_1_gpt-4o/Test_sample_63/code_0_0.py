from itertools import permutations

# Initial grid setup
grid = [
    ['','g','b','','','',''],
    ['','','','','c','',''],
    ['b','','','','','a',''],
    ['f','d','c','','','',''],
    ['','c','e','a','','b',''],
    ['','','','g','','f','d'],
    ['e','a','g','','f','','']
]

# Function to check if a letter can be placed at grid[row][col]
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
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            
            # Try to fill the rest of the grid
            if fill_grid(grid, 0, 0):
                return grid
            
            # Reset the minor diagonal if it doesn't work
            for r in range(7):
                grid[r][6-r] = ''
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