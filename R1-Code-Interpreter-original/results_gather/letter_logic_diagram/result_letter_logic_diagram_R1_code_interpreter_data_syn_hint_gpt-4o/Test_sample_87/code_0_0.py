from itertools import permutations

# Initial grid setup
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

# Function to check if a grid is valid
def is_valid(grid):
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
            return False
    # Check minor diagonal
    minor_diagonal_letter = grid[0][6]
    for i in range(7):
        if grid[i][6-i] != minor_diagonal_letter:
            return False
    return True

# Function to fill the grid
def fill_grid(grid):
    # Try all permutations for the empty rows
    for perm in permutations('abcdefg'):
        # Fill the empty rows
        for i in range(7):
            if '' in grid[i]:
                grid[i] = [perm[j] if cell == '' else cell for j, cell in enumerate(grid[i])]
        # Check if the grid is valid
        if is_valid(grid):
            return grid
    return None

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")