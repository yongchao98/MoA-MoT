def fill_grid(grid):
    # Possible letters
    letters = set('abcdefg')
    
    # Fill the minor diagonal with a letter that can fit
    # Check which letter can be used for the minor diagonal
    for letter in letters:
        if all(grid[i][6-i] == '' or grid[i][6-i] == letter for i in range(7)):
            # Fill the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            break
    
    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find the missing letter for this row and column
                row_letters = set(grid[i])
                col_letters = set(grid[k][j] for k in range(7))
                missing_letter = (letters - row_letters - col_letters).pop()
                grid[i][j] = missing_letter
    
    return grid

# Initial grid
grid = [
    ['g', 'b', 'e', '', 'd', 'f', ''],
    ['b', '', 'c', '', '', 'a', ''],
    ['e', 'c', '', 'f', '', '', 'b'],
    ['', 'd', 'f', '', '', '', ''],
    ['', '', '', '', 'b', 'e', 'c'],
    ['', '', '', '', '', '', 'd'],
    ['', '', 'b', '', '', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))