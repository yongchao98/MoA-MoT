from itertools import permutations

# Initial grid with pre-filled values
grid = [
    ['d', 'f', 'e', 'b', 'a', 'g'],
    ['f', '', 'b', 'a', 'g', ''],
    ['e', 'b', '', '', 'c', 'd'],
    ['b', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
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
        # Place the letter on the minor diagonal
        for i in range(7):
            grid[i][6-i] = letter
        
        # Try to fill the rest of the grid
        if fill_grid(grid, 0, 0):
            return grid

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
for row in solution:
    print(','.join(row))