from itertools import permutations

# Initial grid with empty cells as None
grid = [
    [None, None, None, 'e', 'd', 'c', None],
    ['f', None, 'e', None, None, None, None],
    [None, None, 'd', None, None, 'g', None],
    ['e', None, 'c', None, None, 'b', 'f'],
    ['d', None, None, 'b', None, 'a', 'e'],
    ['c', 'g', 'b', None, 'a', 'e', 'd'],
    [None, 'b', 'f', None, None, 'd', None]
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
    # Determine the letter for the minor diagonal
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Place the letter on the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            break

    # Fill the rest of the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        break

    return grid

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))