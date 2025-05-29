def fill_grid(grid):
    # Determine the letter for the minor diagonal
    # Check which letter is not present in the minor diagonal positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    
    # Find the letter for the minor diagonal
    for letter in letters:
        if all(grid[i][j] != letter for i, j in minor_diagonal_positions):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        row_letters = set(grid[i])
        missing_letters = letters - row_letters
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can be placed in this position
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break
    
    return grid

# Initial grid
grid = [
    ['', '', 'b', '', '', 'd', 'a'],
    ['', '', 'e', '', 'd', '', 'g'],
    ['', 'e', 'c', 'd', '', '', ''],
    ['', 'c', '', '', '', '', ''],
    ['', 'd', 'a', '', 'f', '', 'e'],
    ['', '', 'g', '', '', '', ''],
    ['', '', 'f', '', 'e', 'c', 'd']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))