# Initial grid setup
grid = [
    ['b', '', '', 'g', '', 'c', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', '', 'a', 'c', '', 'b', 'f'],
    ['g', '', '', 'e', '', '', 'd'],
    ['', '', 'e', 'b', '', '', ''],
    ['', 'e', '', 'f', 'd', 'g', 'a'],
    ['e', 'b', 'f', 'd', 'g', '', 'c']
]

# Function to check if a letter can be placed in a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

# Function to fill the grid
def fill_grid(grid):
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = 'a'  # Choose 'a' as it is not over-represented
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    letters = {'a', 'b', 'c', 'd', 'e', 'f', 'g'}
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in letters:
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        break
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))