from itertools import permutations

# Initial grid setup
grid = [
    [' ', 'f', 'b', 'e', 'g', 'd', 'a'],
    [' ', 'b', 'e', 'g', ' ', 'a', 'c'],
    [' ', ' ', ' ', 'd', 'a', ' ', ' '],
    ['e', ' ', ' ', 'a', 'c', 'f', ' '],
    ['g', ' ', 'a', 'c', 'f', 'b', ' '],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', ' ', ' ', 'e', 'g', 'd']
]

# Function to check if placing a letter is valid
def is_valid(grid, row, col, letter):
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
        if all(is_valid(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            
            # Try to fill the rest of the grid
            if fill_grid(grid):
                return grid
            
            # Reset the minor diagonal if it doesn't work
            for i in range(7):
                grid[i][6-i] = ' '
    return None

# Function to fill the grid
def fill_grid(grid):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == ' ':
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if fill_grid(grid):
                            return True
                        grid[row][col] = ' '
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