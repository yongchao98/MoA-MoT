from itertools import permutations

# Initial grid with empty cells as None
grid = [
    [None, None, 'a', None, 'c', 'f', 'b'],
    ['g', 'a', None, None, None, 'b', 'e'],
    [None, None, 'c', 'f', 'b', 'e', 'g'],
    ['d', 'c', 'f', None, 'e', 'g', 'a'],
    [None, 'f', 'b', 'e', 'g', None, None],
    ['f', 'b', 'e', 'g', 'a', 'd', None],
    ['b', 'e', None, 'a', 'd', None, 'f']
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

# Function to fill the grid
def fill_grid(grid):
    # Try to fill the minor diagonal with a single letter
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
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

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))