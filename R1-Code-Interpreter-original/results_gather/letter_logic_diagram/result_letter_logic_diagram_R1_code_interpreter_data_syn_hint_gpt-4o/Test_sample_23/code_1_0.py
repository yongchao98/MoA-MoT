def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve_grid(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve_grid(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != ' ':
        return solve_grid(grid, minor_diagonal_letter, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ' '
    
    return False

def fill_grid(grid):
    for minor_diagonal_letter in 'abcdefg':
        # Set the minor diagonal to the same letter
        for i in range(7):
            grid[i][6-i] = minor_diagonal_letter
        
        # Try to solve the grid with the current minor diagonal letter
        if solve_grid(grid, minor_diagonal_letter):
            return grid

# Initial grid setup
grid = [
    [' ', 'c', 'g', 'b', 'a', 'e', 'f'],
    [' ', ' ', ' ', 'a', 'e', 'f', 'd'],
    ['g', 'b', ' ', 'e', ' ', ' ', 'c'],
    ['b', 'a', 'e', ' ', 'd', 'c', ' '],
    ['a', 'e', 'f', 'd', 'c', ' ', 'b'],
    ['e', 'f', 'd', ' ', 'g', 'b', ' '],
    ['f', ' ', ' ', 'g', ' ', 'a', 'e']
]

# Solve the grid
solution = fill_grid(grid)

# Print the solution
for row in solution:
    print(','.join(row))