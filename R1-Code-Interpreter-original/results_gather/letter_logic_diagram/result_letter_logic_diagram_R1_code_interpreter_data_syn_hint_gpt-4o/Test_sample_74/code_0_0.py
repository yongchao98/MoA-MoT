# Initial grid with empty cells represented by None
grid = [
    [None, None, None, None, 'd', 'g', 'c'],
    ['a', None, None, 'd', None, 'c', None],
    ['a', None, 'd', 'g', None, None, None],
    [None, None, 'g', 'c', None, None, 'a'],
    ['d', None, 'c', 'e', None, None, 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', None, 'd', 'g']
]

# Function to fill the grid
def fill_grid(grid):
    # Determine the letter for the minor diagonal
    # We will try each letter from 'a' to 'g' and see if it fits
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(grid[i][6-i] in (None, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] is None:
                # Find the missing letter for this row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (set('abcdefg') - row_letters - col_letters).pop()
                grid[i][j] = missing_letter

    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))