from itertools import permutations

# Initial grid with empty cells represented as ''
grid = [
    ['f', 'e', 'g', 'a', '', '', ''],
    ['e', 'g', 'a', '', '', '', 'f'],
    ['g', '', '', '', 'b', 'f', 'e'],
    ['', 'd', '', '', 'f', 'e', 'g'],
    ['d', '', 'b', '', 'e', '', 'a'],
    ['', 'b', 'f', 'e', '', 'a', 'd'],
    ['', '', 'e', '', 'a', 'd', 'c']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

# Function to solve the grid
def solve(grid):
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            
            # Try to fill the rest of the grid
            if fill_grid(grid):
                return grid
            
            # Reset the minor diagonal
            for i in range(7):
                grid[i][6-i] = ''
    return None

# Function to fill the grid
def fill_grid(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        if fill_grid(grid):
                            return True
                        grid[row][col] = ''
                return False
    return True

# Solve the grid
solution = solve(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")